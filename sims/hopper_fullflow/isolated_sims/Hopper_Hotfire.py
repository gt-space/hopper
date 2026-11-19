"""Test-stand hotfire: the engine fed from constant-pressure run tanks.

python isolated_sims/Hopper_Hotfire.py [--config case.json] [--out folder] [--check-inputs]

Inputs are in isolated_sims/hotfire_config.json: reservoir set pressures, test-stand
CdAs, throttle-valve full-open CdA, stand-to-injector lines, injectors, regen, the
engine (same fields as the engine section of hopper_config.json) and the thrust
profile. Fluids, the chamber map, ambient pressure and the main-valve sequence
timing come from hopper_config.json.

Network (each side)
    LOX tank  -> LOX stand CdA  -> OTV -> LOX line  ----------------> OINJ -> Main Chamber -> nozzle
    Fuel tank -> Fuel stand CdA -> FTV -> fuel line -> Regen Jacket -> FINJ -> Main Chamber

OTV/FTV are flow controllers (the vehicle ThrottleValve): each opens to whatever CdA
passes the flow the thrust profile needs, up to its full-open CdA. A profile the feed
cannot deliver shows up as a fully open valve and a thrust shortfall.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_CONFIG = ROOT / "isolated_sims" / "hotfire_config.json"
VEHICLE_CONFIG = ROOT / "hopper_config.json"
PSI = 6894.757293168361
LBF = 4.4482216152605
G0 = 9.80665

ASSUMPTIONS = [
    "Run tanks are constant-pressure reservoirs (no blowdown, no pressurant, no ullage).",
    "Each side: tank -> test-stand CdA -> throttle valve -> line -> (fuel: regen) -> injector -> chamber.",
    "The test-stand CdA is one fixed incompressible restriction for all stand plumbing upstream of the throttle "
    "valve; it has a 1 cm passage length for flow inertia only (no steady dP).",
    "OTV/FTV are flow controllers that track the thrust-targeter flows, capped at the full-open CdA.",
    "Engine: chamber map from hopper_config.json; throat, exit area, MR and c*/Cf efficiencies from the hotfire "
    "engine section. Propellant ignites as soon as it reaches the chamber, which has no volume.",
    "Liquids are incompressible at the reservoir pressure and temperature; the regen heats the fuel by the "
    "configured dT before the fuel injector.",
    "Main-valve opening and closing timing (sequence section) come from hopper_config.json.",
    "Not modelled: torch igniter, purges, line priming and chill-down, stand instrumentation.",
]


# ---- Inputs ---- #
def _unique(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON field: {key}")
        value[key] = item
    return value


def _keys(section, name, required, optional=()):
    missing = set(required) - section.keys()
    extra = section.keys() - set(required) - set(optional) - {"note"}
    if missing or extra:
        raise ValueError(f"{name}: missing {sorted(missing)}; unknown {sorted(extra)}.")


def load_hotfire(path):
    path = Path(path)
    config = json.loads(path.read_text(encoding="utf-8-sig"), object_pairs_hook=_unique)
    _keys(config, "hotfire config", ("simulation", "engine", "profile", "oxidizer", "fuel"), ("description",))
    _keys(config["simulation"], "simulation", ("maximum_step_s", "save_interval_s", "output_directory"))
    _keys(config["profile"], "profile", ("csv", "interpolation"))
    side_keys = ("reservoir_pressure_psia", "temperature_k", "stand_cda_mm2", "throttle_valve_max_cda_mm2",
                 "throttle_valve_actuation_s", "line", "injector_cda_mm2", "injector_length_in")
    _keys(config["oxidizer"], "oxidizer", side_keys)
    _keys(config["fuel"], "fuel", side_keys + ("regen",))
    for side in ("oxidizer", "fuel"):
        spec = config[side]
        _keys(spec["line"], f"{side}.line", ("length_in", "outer_diameter_in", "wall_thickness_in", "roughness_in"))
        for key in ("reservoir_pressure_psia", "temperature_k", "stand_cda_mm2", "throttle_valve_max_cda_mm2",
                    "injector_cda_mm2"):
            if spec[key] <= 0:
                raise ValueError(f"{side}.{key} must be positive.")
        line = spec["line"]
        if line["length_in"] <= 0 or not 0 < 2 * line["wall_thickness_in"] < line["outer_diameter_in"]:
            raise ValueError(f"{side}.line: need a positive length and a wall thinner than half the OD.")
    _keys(config["fuel"]["regen"], "fuel.regen", ("cda_mm2", "pc_psia_delta_t_k", "length_in"))
    profile = config["profile"]
    if profile["interpolation"] not in ("linear", "step"):
        raise ValueError("profile.interpolation must be linear or step.")
    csv_path = (path.parent / profile["csv"]).resolve()
    with csv_path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["time [s]", "thrust [lbf]"]:
            raise ValueError(f"{csv_path.name}: columns must be 'time [s]', 'thrust [lbf]'.")
        points = [(float(row["time [s]"]), float(row["thrust [lbf]"])) for row in reader]
    if len(points) < 2 or points[0][0] != 0 or any(b[0] <= a[0] for a, b in zip(points, points[1:])) \
            or any(f < 0 for _, f in points) or not any(f > 0 for _, f in points):
        raise ValueError(f"{csv_path.name}: start at t = 0, times increasing, thrusts >= 0 and not all zero.")
    profile["csv_path"] = str(csv_path)   # in memory only, for the run
    profile["points_s_lbf"] = points
    return config


def thrust_targets(points):
    """Targeter Pc and flows for every profile point. A 0 lbf point has both
    throttle valves shut: zero flow targets, Pc at ambient. The targeter needs the
    first point it solves at t = 0, so the positive points are shifted to start
    there and shifted back."""
    import HOPPER_DATA as D
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
    positive = [(t, f) for t, f in points if f > 0]
    shift = positive[0][0]
    targets = HopperThrustTargets({t - shift: f * LBF for t, f in positive})
    for target in targets:
        target["time"] += shift
    targets += [{"time": t, "thrust": 0.0, "chamber_pressure": D.AMBIENT_PRESSURE, "ox_mass_flow": 0.0,
                 "fuel_mass_flow": 0.0} for t, f in points if f == 0]
    return sorted(targets, key=lambda target: target["time"])


def feed_elements(config):
    """Element specs (hopper_config.json component types) for the two feeds."""
    def side(spec, label, valve, injector, line_outlet, regen=None):
        elements = [
            {"type": "Orifice", "name": f"{label} Stand", "cda_mm2": spec["stand_cda_mm2"],
             "length_in": 0.3937007874015748, "outlet_node": f"{valve} Inlet"},
            {"type": "ThrottleValve", "name": valve, "max_cda_mm2": spec["throttle_valve_max_cda_mm2"],
             "outlet_node": f"{valve} Outlet", "open_time_s": spec["throttle_valve_actuation_s"],
             "close_time_s": spec["throttle_valve_actuation_s"]},
            {"type": "Line", "name": f"{label} Line", "height_change_in": 0.0, "outlet_node": line_outlet,
             **spec["line"]},
        ]
        if regen:
            elements.append({"type": "Regen", "name": "Regen Jacket", "outlet_node": "Fuel Injector Manifold", **regen})
        elements.append({"type": "Injector", "name": injector, "cda_mm2": spec["injector_cda_mm2"],
                         "length_in": spec["injector_length_in"], "outlet_node": f"{injector} Outlet"})
        return elements

    return {"LOX": side(config["oxidizer"], "LOX", "OTV", "OINJ", "LOX Injector Manifold"),
            "Fuel": side(config["fuel"], "Fuel", "FTV", "FINJ", "Regen Inlet Manifold", config["fuel"]["regen"])}


def vehicle_config_with_engine(config, config_dir, output):
    """hopper_config.json with its engine section replaced by the hotfire engine,
    written to the output folder. Relative paths are made absolute so the copy
    loads from there. HOPPER_DATA (and everything built on it) reads this copy."""
    vehicle = json.loads(VEHICLE_CONFIG.read_text(encoding="utf-8-sig"))
    engine = copy.deepcopy(config["engine"])
    engine.pop("note", None)
    for key in ("cstar_efficiency", "cf_efficiency"):
        if engine[key].get("csv"):
            engine[key]["csv"] = str((config_dir / engine[key]["csv"]).resolve())
    vehicle["engine"] = engine
    for key, value in vehicle["files"].items():
        vehicle["files"][key] = str((ROOT / value).resolve())
    for spec in vehicle["thrust_profile"]["profiles"].values():
        if spec.get("csv"):
            spec["csv"] = str((ROOT / spec["csv"]).resolve())
    vehicle["monte_carlo"]["directory"] = str((ROOT / vehicle["monte_carlo"]["directory"]).resolve())
    path = output / "hopper_config_hotfire_engine.json"
    path.write_text(json.dumps(vehicle, indent=2) + "\n", encoding="utf-8")
    return path


# ---- Model ---- #
def run(config, output, check_inputs=False):
    import HOPPER_DATA as D
    from fullflow import Network, Sequence, State, Volume
    from scipy.constants import psi
    import Backend.Hopper_Prop_System as prop
    from Backend.Hopper_Components import EnginePerformance
    from Backend.Hopper_Config import build_component
    from Backend.Hopper_Control import FunctionLookup
    from Backend.Hopper_Solver import transient_solver

    print("Solving thrust targets...", flush=True)
    targets = thrust_targets(config["profile"]["points_s_lbf"])
    profile_start, profile_end, shutdown_end = prop.burn_timeline(targets)
    t_final = shutdown_end + D.TAIL_OFF_TIME
    timeline = {"main_fuel_opens_s": D.MAIN_FUEL_OPEN_TIME, "main_ox_opens_s": D.MAIN_OX_OPEN_TIME,
                "profile_start_s": profile_start, "profile_end_s": profile_end,
                "shutdown_end_s": shutdown_end, "end_s": t_final}
    hardware = {label: [build_component(copy.deepcopy(spec)) for spec in specs]
                for label, specs in feed_elements(config).items()}
    if check_inputs:
        print(json.dumps(timeline, indent=2))
        for label, elements in hardware.items():
            print(label, " -> ".join(e.name for e in elements), "-> Main Chamber")
        return None, targets, timeline, hardware

    # build_liquid_feed reads the liquid temperature (and, for regen, the property
    # pressure) from the vehicle tank objects: point them at the reservoirs.
    sides = {"LOX": (config["oxidizer"], D.LOX_TANK), "Fuel": (config["fuel"], D.FUEL_TANK)}
    for spec, tank in sides.values():
        tank.liquid_temperature = spec["temperature_k"]
        tank.set_pressure = spec["reservoir_pressure_psia"] * psi

    start = time.perf_counter()
    print("Building the test-stand network...", flush=True)
    net = Network("Hopper Hotfire")
    flow_targets = {"LOX": State(0.0), "Fuel": State(0.0)}
    thrust_target = State(0.0)
    interpolation = config["profile"]["interpolation"]
    for target, key, name in [(flow_targets["LOX"], "ox_mass_flow", "OTV Flow Target"),
                              (flow_targets["Fuel"], "fuel_mass_flow", "FTV Flow Target"),
                              (thrust_target, "thrust", "Thrust Target")]:
        times, values = prop.target_schedule(targets, key, interpolation)
        Sequence(name, net, target=target, times=times, values=values)

    PC = State(D.AMBIENT_PRESSURE)
    feeds = {}
    for label, (spec, tank) in sides.items():
        pressure = spec["reservoir_pressure_psia"] * psi
        liquid = FunctionLookup(f"{label} Tank Liquid", net, prop.liquid_properties, propellant=tank.propellant,
                                pressure=pressure, temperature=spec["temperature_k"], outputs=prop.LIQUID_PROPERTIES)
        reservoir = {"liquid": liquid, "liquid_node": SimpleNamespace(pressure=State(pressure))}
        feeds[label] = prop.build_liquid_feed(net, label, hardware[label], reservoir, PC, flow_targets[label])
        net.track(f"{label} Tank Pressure [Pa]", reservoir["liquid_node"].pressure)

    floor = D.MIXTURE_RATIO_FLOW_FLOOR
    ox_in, fuel_in = feeds["LOX"]["chamber_inflow"], feeds["Fuel"]["chamber_inflow"]
    mixture_ratio = ((ox_in + D.MIXTURE_RATIO * floor) / (fuel_in + floor)).clip(
        D.CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT, D.CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT)
    engine = EnginePerformance("Main Chamber", net, chamber_pressure=PC, mixture_ratio=mixture_ratio,
                               sampled_efficiency=True, chamber_map_filename=D.CHAMBER_MAP_FILENAME)
    chamber = Volume("Main Chamber", net, pressure=PC, mass_flow_in=ox_in + fuel_in)
    chamber.mass_flow_out = engine["mass_flow"]

    net.track("Chamber Pressure [psia]", PC / psi)
    net.track("Thrust [lbf]", engine["thrust"] / LBF)
    net.track("Target Thrust [lbf]", thrust_target / LBF)
    net.track("Mixture Ratio", mixture_ratio)
    net.track("Chamber Temperature [K]", engine["gas"].temperature)
    net.track("Nozzle Mass Flow [kg/s]", engine["mass_flow"])
    net.track("C* Efficiency", engine["cstar_efficiency"])
    net.track("Cf Efficiency", engine["cf_efficiency"])
    net.track("Nozzle Exit Pressure [psia]", engine["nozzle"].exit_static_pressure / psi)

    dt = config["simulation"]["maximum_step_s"]
    filename = output / "Hopper_Hotfire.h5"
    print(f"Network ready in {time.perf_counter() - start:.1f} s; solving {t_final:g} s at <= {dt:g} s steps.",
          flush=True)
    transient_solver(net, max_evaluations=D.PROP_MAX_NFEV).solve(
        dt=dt, t_final=t_final, save_dt=config["simulation"]["save_interval_s"], filename=str(filename),
        rtol=D.PROP_RTOL, jacobian_method=D.PROP_JACOBIAN_METHOD, verbose=False, statistics=True)
    print(f"Solved in {time.perf_counter() - start:.1f} s.", flush=True)
    return filename, targets, timeline, hardware


# ---- Results ---- #
def read_history(filename, hardware):
    """Time histories in engineering units: chamber, then each element of each feed
    in flow order (pressures, dP, flow), throttle-valve CdA/opening, injector
    stiffness and regen heating."""
    import h5py
    import numpy as np
    from Backend.Hopper_Plots import _safe_name

    with h5py.File(filename, "r") as h5:
        network = next(name for name, group in h5.items() if isinstance(group, h5py.Group) and "transient" in group)
        base = h5[f"{network}/transient/runs/base"]
        tracks = base["tracks"]

        def track(name, scale=1.0):
            return np.asarray(tracks[_safe_name(name)], dtype=float) * scale

        h = {"Time [s]": np.asarray(base["time"], dtype=float)}
        for name in ("Thrust [lbf]", "Target Thrust [lbf]", "Chamber Pressure [psia]", "Mixture Ratio",
                     "Chamber Temperature [K]", "C* Efficiency", "Cf Efficiency", "Nozzle Exit Pressure [psia]"):
            h[name] = track(name)
        pc = h["Chamber Pressure [psia]"]
        for label, elements in hardware.items():
            h[f"{label} Tank Pressure [psia]"] = track(f"{label} Tank Pressure [Pa]", 1 / PSI)
            h[f"{label} Mass Flow [kg/s]"] = track(f"{elements[-1].name} Mass Flow [kg/s]")
            for e in elements:
                h[f"{e.name} Inlet Pressure [psia]"] = track(f"{e.name} Upstream Pressure [Pa]", 1 / PSI)
                h[f"{e.name} dP [psi]"] = track(f"{e.name} dP [psi]")
            valve, injector = elements[1], elements[-1]
            h[f"{valve.name} Target Mass Flow [kg/s]"] = track(f"{valve.name} Target Mass Flow [kg/s]")
            h[f"{valve.name} CdA [mm2]"] = track(f"{valve.name} CdA [mm2]")
            h[f"{valve.name} Open [%]"] = track(f"{valve.name} Opening", 100.0)
            with np.errstate(divide="ignore", invalid="ignore"):
                h[f"{injector.name} Stiffness [%]"] = np.where(pc > 50, 100 * h[f"{injector.name} dP [psi]"] / pc,
                                                               np.nan)
        h["Regen Jacket dT [K]"] = track("Regen Jacket dT [K]")
        h["Fuel Injector Inlet Temperature [K]"] = track("Regen Jacket Outlet Temperature [K]")
    total = h["LOX Mass Flow [kg/s]"] + h["Fuel Mass Flow [kg/s]"]
    with np.errstate(divide="ignore", invalid="ignore"):
        h["Isp [s]"] = np.where(total > 0.05, h["Thrust [lbf]"] * LBF / (total * G0), np.nan)
    return h


def holds(config, timeline):
    """(start, end, thrust) on the sim clock for each constant-thrust stretch of the profile."""
    points = config["profile"]["points_s_lbf"]
    t0 = timeline["profile_start_s"]
    return [(t0 + a, t0 + b, f) for (a, f), (b, g) in zip(points, points[1:]) if f == g and b > a]


def summarize(h, config, timeline, hardware):
    import numpy as np
    t = h["Time [s]"]
    burning = (t >= timeline["profile_start_s"]) & (t <= timeline["profile_end_s"])

    def snapshot(mask):
        return {k: float(np.nanmean(v[mask])) for k, v in h.items() if k != "Time [s]"}

    hold_rows = []
    for a, b, thrust in holds(config, timeline):
        mask = (t >= a + 0.5 * (b - a)) & (t <= b)   # second half of the hold: settled
        if mask.any():
            hold_rows.append({"hold": f"{thrust:g} lbf, {a:.2f}-{b:.2f} s", **snapshot(mask)})

    shortfall = h["Target Thrust [lbf]"] - h["Thrust [lbf]"]
    saturated = {}
    for valve in ("OTV", "FTV"):
        full = burning & (h[f"{valve} Open [%]"] >= 99.5)
        saturated[valve] = {"fully_open": bool(full.any()),
                            "first_full_open_s": float(t[full][0]) if full.any() else None,
                            "target_thrust_when_first_full_open_lbf":
                                float(h["Target Thrust [lbf]"][full][0]) if full.any() else None}
    peak = int(np.argmax(np.where(burning, h["Thrust [lbf]"], -np.inf)))
    dts = np.diff(t)

    def used(column):
        v = h[column]
        return float(np.sum(0.5 * (v[1:] + v[:-1]) * dts))

    return {
        "timeline_s": timeline,
        "max_thrust_during_profile": {"Thrust [lbf]": float(h["Thrust [lbf]"][peak]), "time_s": float(t[peak]),
                                      **{k: float(h[k][peak]) for k in (
                                          "Chamber Pressure [psia]", "Mixture Ratio", "OTV Open [%]", "FTV Open [%]",
                                          "OTV CdA [mm2]", "FTV CdA [mm2]", "OINJ Stiffness [%]", "FINJ Stiffness [%]",
                                          "LOX Mass Flow [kg/s]", "Fuel Mass Flow [kg/s]")}},
        "largest_thrust_shortfall_during_profile_lbf": float(np.max(shortfall[burning])) if burning.any() else None,
        "throttle_valves_fully_open": saturated,
        "minimum_during_profile": {k: float(np.nanmin(h[k][burning])) for k in (
            "OINJ Stiffness [%]", "FINJ Stiffness [%]", "OTV dP [psi]", "FTV dP [psi]")},
        "maximum_during_profile": {k: float(np.nanmax(h[k][burning])) for k in (
            "Chamber Pressure [psia]", "OTV dP [psi]", "FTV dP [psi]", "OTV Open [%]", "FTV Open [%]",
            "OINJ dP [psi]", "FINJ dP [psi]", "Regen Jacket dT [K]")},
        "holds": hold_rows,
        "propellant_used_kg": {"LOX": used("LOX Mass Flow [kg/s]"), "Fuel": used("Fuel Mass Flow [kg/s]")},
        "feeds": {label: [e.name for e in elements] for label, elements in hardware.items()},
        "assumptions": ASSUMPTIONS,
    }


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(h, config, timeline, hardware, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    t = h["Time [s]"]
    burn = (t >= timeline["profile_start_s"]) & (t <= timeline["profile_end_s"])
    ox, fu = "C0", "C3"

    def finish(ax, ylabel, legend=True):
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        if legend:
            ax.legend(fontsize=8)

    # Time histories.
    panels = [
        ("Thrust [lbf]", [("Thrust [lbf]", "Thrust", "k", "-"), ("Target Thrust [lbf]", "Target", "0.5", ":")]),
        ("Chamber pressure [psia]", [("Chamber Pressure [psia]", "Pc", "k", "-")]),
        ("Mass flow [kg/s]", [("LOX Mass Flow [kg/s]", "LOX", ox, "-"), ("OTV Target Mass Flow [kg/s]", "LOX target", ox, ":"),
                              ("Fuel Mass Flow [kg/s]", "Fuel", fu, "-"), ("FTV Target Mass Flow [kg/s]", "Fuel target", fu, ":")]),
        ("Mixture ratio", [("Mixture Ratio", "MR", "k", "-")]),
        ("Throttle valve open [%]", [("OTV Open [%]", "OTV", ox, "-"), ("FTV Open [%]", "FTV", fu, "-")]),
        ("Throttle valve CdA [mm²]", [("OTV CdA [mm2]", "OTV", ox, "-"), ("FTV CdA [mm2]", "FTV", fu, "-")]),
        ("Throttle valve dP [psi]", [("OTV dP [psi]", "OTV", ox, "-"), ("FTV dP [psi]", "FTV", fu, "-")]),
        ("Injector dP [psi]", [("OINJ dP [psi]", "OINJ", ox, "-"), ("FINJ dP [psi]", "FINJ", fu, "-")]),
        ("Injector stiffness [% of Pc]", [("OINJ Stiffness [%]", "OINJ", ox, "-"), ("FINJ Stiffness [%]", "FINJ", fu, "-")]),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 2.3 * len(panels)), sharex=True)
    for ax, (ylabel, lines) in zip(axes, panels):
        for column, name, color, style in lines:
            ax.plot(t, h[column], color=color, ls=style, lw=1.5, label=name)
        finish(ax, ylabel, legend=len(lines) > 1)
    axes[0].set_title("Test-stand hotfire")
    axes[-1].set_xlabel("Time [s]")
    fig.tight_layout()
    fig.savefig(output / "hotfire_overview.png", dpi=140)
    plt.close(fig)

    # Steady maps against thrust (profile window only).
    x = h["Thrust [lbf]"][burn]
    order = np.argsort(x)
    maps = [
        ("Chamber pressure [psia]", [("Chamber Pressure [psia]", "Pc", "k")]),
        ("Throttle valve open [%]", [("OTV Open [%]", "OTV", ox), ("FTV Open [%]", "FTV", fu)]),
        ("Throttle valve dP [psi]", [("OTV dP [psi]", "OTV", ox), ("FTV dP [psi]", "FTV", fu)]),
        ("Injector stiffness [% of Pc]", [("OINJ Stiffness [%]", "OINJ", ox), ("FINJ Stiffness [%]", "FINJ", fu)]),
        ("Mass flow [kg/s]", [("LOX Mass Flow [kg/s]", "LOX", ox), ("Fuel Mass Flow [kg/s]", "Fuel", fu)]),
        ("Isp [s]", [("Isp [s]", "Isp", "k")]),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    for ax, (ylabel, lines) in zip(axes.flat, maps):
        for column, name, color in lines:
            ax.plot(x[order], h[column][burn][order], color=color, lw=1.6, label=name)
        finish(ax, ylabel, legend=len(lines) > 1)
    for ax in axes[-1]:
        ax.set_xlabel("Thrust [lbf]")
    fig.suptitle("Test-stand hotfire: steady values against thrust (profile window)")
    fig.tight_layout()
    fig.savefig(output / "hotfire_vs_thrust.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x[order], h["Chamber Pressure [psia]"][burn][order], color="k", lw=1.8)
    finish(ax, "Chamber pressure [psia]", legend=False)
    ax.set_xlabel("Thrust [lbf]")
    ax.set_title("Chamber pressure vs thrust")
    fig.tight_layout()
    fig.savefig(output / "hotfire_pc_vs_thrust.png", dpi=140)
    plt.close(fig)

    # Element dP and pressure ladder at the peak thrust.
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, (label, elements) in zip(axes, hardware.items()):
        for e in elements:
            ax.plot(t, h[f"{e.name} dP [psi]"], lw=1.5, label=e.name)
        finish(ax, f"{label} element dP [psi]")
    axes[-1].set_xlabel("Time [s]")
    fig.tight_layout()
    fig.savefig(output / "hotfire_element_dp.png", dpi=140)
    plt.close(fig)

    # Stations share x positions across the two feeds; the LOX side skips the regen slot
    # so both end at the chamber.
    peak = int(np.argmax(np.where(burn, h["Thrust [lbf]"], -np.inf)))
    stations = ["Tank", "Stand out", "Valve out", "Line out", "Regen out", "Chamber"]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    for (label, elements), color, offset in zip(hardware.items(), (ox, fu), (-16, 8)):
        x = [0, 1, 2, 3] + ([4] if len(elements) == 5 else []) + [5]
        pressures = [h[f"{label} Tank Pressure [psia]"][peak]] + \
                    [h[f"{e.name} Inlet Pressure [psia]"][peak] - h[f"{e.name} dP [psi]"][peak] for e in elements]
        ax.plot(x, pressures, "o-", color=color, label=label)
        for i, p in zip(x, pressures):
            ax.annotate(f"{p:.0f}", (i, p), textcoords="offset points", xytext=(0, offset), ha="center",
                        fontsize=8, color=color)
    ax.set_xticks(range(len(stations)), stations)
    finish(ax, "Pressure [psia]")
    ax.set_title(f"Pressure ladder at peak thrust ({h['Thrust [lbf]'][peak]:.0f} lbf, t = {t[peak]:.2f} s)")
    fig.tight_layout()
    fig.savefig(output / "hotfire_pressure_ladder.png", dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Test-stand hotfire from constant-pressure run tanks.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG), help="hotfire JSON (default isolated_sims/hotfire_config.json)")
    parser.add_argument("--out", help="output folder (default <simulation.output_directory>/<UTC time>)")
    parser.add_argument("--check-inputs", action="store_true", help="build the hardware and print the timeline only")
    args = parser.parse_args()
    started = time.perf_counter()
    config_path = Path(args.config).resolve()
    config = load_hotfire(config_path)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    output = Path(args.out) if args.out else ROOT / config["simulation"]["output_directory"] / stamp
    if args.check_inputs and not args.out:
        import tempfile
        output = Path(tempfile.mkdtemp(prefix="hotfire_check_"))
    output.mkdir(parents=True, exist_ok=bool(args.check_inputs))
    # HOPPER_DATA reads HOPPER_CONFIG when first imported (inside run), so the
    # engine everything is built with is the hotfire engine.
    os.environ["HOPPER_CONFIG"] = str(vehicle_config_with_engine(config, config_path.parent, output))

    if args.check_inputs:
        run(config, output, check_inputs=True)
        return
    # Copy the inputs next to the results; the copied config points at the copied profile.
    profile_csv = Path(config["profile"]["csv_path"])
    shutil.copy2(profile_csv, output / profile_csv.name)
    snapshot = json.loads(config_path.read_text(encoding="utf-8-sig"))
    snapshot["profile"]["csv"] = profile_csv.name
    (output / "hotfire_config.json").write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
    filename, targets, timeline, hardware = run(config, output)
    h = read_history(filename, hardware)
    write_csv(output / "hotfire_history.csv", [dict(zip(h, row)) for row in zip(*h.values())])
    write_csv(output / "thrust_targets.csv", [
        {"profile time [s]": p["time"], "thrust [lbf]": p["thrust"] / LBF,
         "Pc [psia]": p["chamber_pressure"] / PSI, "LOX [kg/s]": p["ox_mass_flow"], "fuel [kg/s]": p["fuel_mass_flow"]}
        for p in targets])
    summary = summarize(h, config, timeline, hardware)
    if summary["holds"]:
        write_csv(output / "hotfire_holds.csv", summary["holds"])
    summary["elapsed_s"] = time.perf_counter() - started
    (output / "hotfire_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    plot(h, config, timeline, hardware, output)
    print(json.dumps({k: summary[k] for k in ("max_thrust_during_profile", "throttle_valves_fully_open",
                                              "minimum_during_profile", "maximum_during_profile")}, indent=2))
    print(f"Results in {output}")


if __name__ == "__main__":
    main()

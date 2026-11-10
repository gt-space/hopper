"""Test-stand hotfire: the engine fed from constant-pressure run tanks.

python isolated_sims/Hopper_Hotfire.py [--config case.json] [--out folder] [--check-inputs]

Inputs are in isolated_sims/hotfire_config.json (hardware, reservoirs, thrust ramp,
torch and purge timing). The engine (throat, nozzle, c*/Cf efficiencies, chamber
map, target MR), the fluids and the OTV/FTV opening timing come from
hopper_config.json, so the engine is the same one the vehicle sim uses.

Network (test-stand P&ID)
    LOX reservoir      -> LOX line -> OTV -> line -> OINJ ------------------> Main Chamber -> nozzle
    Kerosene reservoir -> fuel line -> FTV -> line -> regen -> FINJ --------> Main Chamber
    OTV Inlet -> OIMV -> OICK -> OIORF -> OTHX -> torch ox injector  --> Torch Chamber
    FTV Inlet -> FIMV -> FICK -> FIORF -> torch fuel injector       --> Torch Chamber -> torch throat -> Main Chamber
    GN2 supply -> Purge Reg -> Purge RV (relief to ambient)
                            -> OPV  -> OPCK  -> below OTV
                            -> FPV  -> FPCK  -> below FTV
                            -> FIPV -> FIPCK -> FIPORF -> torch fuel injector -> Torch Chamber

OTV/FTV are flow controllers that follow the thrust targets (same ThrottleValve
as the vehicle). The torch chamber is a zero-volume combustion node on the same
chamber map; its products flow through the torch throat (ideal, bidirectional,
subsonic or choked) into the main chamber and count in the main chamber's flow and
MR. Purge gas only sees its back pressure and leaves the model there (as in the
vehicle sim). Assumptions are listed in ASSUMPTIONS below and in the run summary.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_CONFIG = ROOT / "isolated_sims" / "hotfire_config.json"

ASSUMPTIONS = [
    "Run tanks are constant-pressure reservoirs (no blowdown, no pressurant, no ullage).",
    "Engine, chamber map, efficiencies, target MR and the OTV/FTV opening ramp come from hopper_config.json.",
    "OTV/FTV are flow controllers that track the thrust-targeter flows; the targets ignore the torch flow.",
    "Propellant ignites as soon as it reaches a chamber; both chambers have no volume.",
    "Torch: same LOX/RP-1 chamber map, c* efficiency from hotfire_config.json, ideal bidirectional throat.",
    "OTHX is a plain tube (friction-equivalent CdA): no heating or vaporization of the LOX.",
    "Liquids are incompressible at the reservoir pressure and temperature; regen heating from hopper_config.json.",
    "Purge gas: ideal regulator off a constant GN2 supply; the gas only sees its back pressure and is not added "
    "to the propellant lines or chambers. Check valves block reverse flow.",
    "Torch purge path includes the torch fuel injector; it is solved as if FIMV is shut (no shared-node coupling).",
]


def _unique(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON field: {key}")
        value[key] = item
    return value


def load_hotfire(path):
    path = Path(path)
    config = json.loads(path.read_text(encoding="utf-8-sig"), object_pairs_hook=_unique)
    required = {"simulation", "reservoirs", "thrust_profile", "torch", "purges", "purge_supply", "purge_regulator",
                "oxidizer_feed", "fuel_feed", "torch_oxidizer_leg", "torch_fuel_leg", "purge_branches"}
    missing, extra = required - config.keys(), config.keys() - required - {"description", "notes"}
    if missing or extra:
        raise ValueError(f"{path.name}: missing sections {sorted(missing)}; unknown sections {sorted(extra)}.")
    torch, purges = config["torch"], config["purges"]
    if not 0 <= torch["open_s"] < torch["close_s"] or torch["valve_stroke_s"] <= 0 \
            or torch["close_s"] - torch["open_s"] < 2 * torch["valve_stroke_s"]:
        raise ValueError("torch: need 0 <= open_s < close_s with room for two valve strokes.")
    if not 0 < torch["cstar_efficiency"] <= 1 or torch["throat_area_in2"] <= 0:
        raise ValueError("torch: throat area must be positive and c* efficiency in (0, 1].")
    if purges["main_duration_s"] <= 0 or purges["valve_stroke_s"] <= 0:
        raise ValueError("purges: duration and valve stroke must be positive.")
    ramp = config["thrust_profile"]
    if ramp["interpolation"] not in ("linear", "step") or ramp["hold_s"] < 0:
        raise ValueError("thrust_profile: interpolation linear/step and a nonnegative hold.")
    for side in ("oxidizer", "fuel"):
        if config["reservoirs"][side]["pressure_psia"] <= 0:
            raise ValueError(f"reservoirs.{side}: pressure must be positive.")
    return config


def window(open_time, close_time, stroke):
    """Valve schedule (times, openings): closed, ramp open, hold, ramp closed."""
    return ([0.0, open_time, open_time + stroke, close_time, close_time + stroke], [0.0, 0.0, 1.0, 1.0, 0.0])


def thrust_profile(spec, lbf):
    from Backend.Hopper_Inputs import ramp_profile
    profile = ramp_profile(spec["start_lbf"] * lbf, spec["end_lbf"] * lbf, spec["ramp_duration_s"], spec["points"])
    if spec["hold_s"] > 0:
        profile[spec["ramp_duration_s"] + spec["hold_s"]] = spec["end_lbf"] * lbf
    return profile


def build_hardware(config, schedules):
    """Element objects from the JSON specs (same component types as hopper_config.json)."""
    from Backend.Hopper_Config import build_component, build_branch

    def component(spec):
        spec = copy.deepcopy(spec)
        if spec["type"] == "Tube":   # liquid tube as its friction-equivalent CdA, A / sqrt(1 + f L / D)
            spec["type"] = "GasLine"
        return build_component(spec, schedules)

    hardware = {"oxidizer_feed": [component(s) for s in config["oxidizer_feed"]],
                "fuel_feed": [component(s) for s in config["fuel_feed"]],
                "regulator": component(config["purge_regulator"]),
                "branches": [build_branch(copy.deepcopy(s), schedules) for s in config["purge_branches"]]}
    for key in ("torch_oxidizer_leg", "torch_fuel_leg"):
        leg = config[key]
        hardware[key] = {"source": leg["source"], "elements": [component(s) for s in leg["elements"]]}
    names = [e.name for k in ("oxidizer_feed", "fuel_feed") for e in hardware[k]]
    names += [e.name for k in ("torch_oxidizer_leg", "torch_fuel_leg") for e in hardware[k]["elements"]]
    names += [e.name for b in hardware["branches"] for e in b.elements] + [b.name for b in hardware["branches"]]
    if len(set(names)) != len(names):
        raise ValueError(f"Component names must be unique: {sorted({n for n in names if names.count(n) > 1})}")
    return hardware


def torch_throat_flow(p_up, p_down, gamma_up, r_up, t_up, gamma_down, r_down, t_down, area, band):
    """Ideal isentropic flow [kg/s] through the torch throat, positive torch -> main.
    Subsonic or choked; the upstream side sets the gas. Linear in dP inside +/- band
    so the zero-flow state (torch off) is smooth."""
    def one_way(p0, pb, g, r, t):
        ratio = max(pb / p0, (2.0 / (g + 1.0)) ** (g / (g - 1.0)))
        term = ratio ** (2.0 / g) - ratio ** ((g + 1.0) / g)
        return area * p0 * math.sqrt(max(term, 0.0) * 2.0 * g / ((g - 1.0) * r * t))

    dp = p_up - p_down
    if abs(dp) < band:
        edge = one_way(max(p_up, p_down) + band - abs(dp), min(p_up, p_down), gamma_up, r_up, t_up)
        return {"mass_flow": edge * dp / band}
    if dp > 0:
        return {"mass_flow": one_way(p_up, p_down, gamma_up, r_up, t_up)}
    return {"mass_flow": -one_way(p_down, p_up, gamma_down, r_down, t_down)}


def build_torch_leg(network, name, elements, source_pressure, fluid, sink_pressure, track_element):
    """Liquid restrictions in series (SeriesCdA) from a feed node into the torch
    chamber; scheduled valves scale the discharge coefficient; no reverse flow."""
    from fullflow import DischargeCoefficient, SeriesCdA, State
    from Backend.Hopper_Elements import FeedContext

    context = FeedContext(network, fluid, name)
    openings = [element.opening(context) for element in elements]
    discharge = 1.0
    for opening in openings:
        discharge = discharge * opening
    path = SeriesCdA(f"{name} Series CdA", network, effective_areas=[e.cda_value for e in elements])
    component = DischargeCoefficient(name, network, upstream_pressure=source_pressure, downstream_pressure=sink_pressure,
                                     density=fluid.density, discharge_coefficient=discharge,
                                     cross_sectional_area=path.effective_area,
                                     mass_flow=State(0.0))   # assigned before the first solve (MR reads it)
    mass_flow = State.maximum(component.mass_flow, 0.0)
    network.track(f"{name} Mass Flow [kg/s]", mass_flow)
    network.track(f"{name} Effective CdA [m2]", path.effective_area * discharge)
    for element, opening in zip(elements, openings):
        element._connect(source_pressure, sink_pressure)
        element.mass_flow, element.cda, element.state, element.group = mass_flow, element.cda_value, opening, name
        track_element(network, element, fluid.enthalpy)
    return mass_flow


def run(config, output, check_inputs=False):
    import HOPPER_DATA as D
    from fullflow import Lookup, Map, Network, Sequence, State, Volume
    from thermoprop import Fluid
    from scipy.constants import psi, milli, pound_force
    import Backend.Hopper_Prop_System as prop
    from Backend.Hopper_Components import EnginePerformance
    from Backend.Hopper_Control import FunctionLookup
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
    from Backend.Hopper_Solver import transient_solver

    # ---- Thrust targets and timeline ---- #
    print("Solving thrust targets...", flush=True)
    targets = HopperThrustTargets(thrust_profile(config["thrust_profile"], pound_force))
    profile_start, profile_end, shutdown_end = prop.burn_timeline(targets)
    torch, purges = config["torch"], config["purges"]
    purge_open = shutdown_end + purges["main_open_after_shutdown_s"]
    purge_close = purge_open + purges["main_duration_s"]
    torch_purge_close = shutdown_end + purges["torch_purge_close_after_shutdown_s"]
    if torch_purge_close <= torch["close_s"] + purges["valve_stroke_s"]:
        raise ValueError("purges.torch_purge_close_after_shutdown_s closes FIPV before it has opened.")
    schedules = {"torch": window(torch["open_s"], torch["close_s"], torch["valve_stroke_s"]),
                 "torch_purge": window(torch["close_s"], torch_purge_close, purges["valve_stroke_s"]),
                 "main_purge": window(purge_open, purge_close, purges["valve_stroke_s"])}
    t_final = max(shutdown_end, purge_close, torch_purge_close) + purges["valve_stroke_s"] \
        + config["simulation"]["tail_off_s"]
    timeline = {"torch_open_s": torch["open_s"], "torch_close_s": torch["close_s"],
                "main_fuel_opens_s": D.MAIN_FUEL_OPEN_TIME, "main_ox_opens_s": D.MAIN_OX_OPEN_TIME,
                "profile_start_s": profile_start,
                "profile_end_s": profile_end, "shutdown_end_s": shutdown_end,
                "main_purges_open_s": purge_open, "main_purges_close_s": purge_close,
                "torch_purge_open_s": torch["close_s"], "torch_purge_close_s": torch_purge_close, "end_s": t_final}
    hardware = build_hardware(config, schedules)
    if check_inputs:
        print(json.dumps(timeline, indent=2))
        return None, targets, timeline

    # build_liquid_feed reads the liquid temperature (and, for regen, the property
    # pressure) from the vehicle tank objects: point them at the reservoirs.
    reservoirs = config["reservoirs"]
    for tank, side in ((D.LOX_TANK, "oxidizer"), (D.FUEL_TANK, "fuel")):
        tank.liquid_temperature = reservoirs[side]["temperature_k"]
        tank.set_pressure = reservoirs[side]["pressure_psia"] * psi

    start = time.perf_counter()
    print("Building the test-stand network...", flush=True)
    net = Network("Hopper Hotfire")
    lox_target, fuel_target, thrust_target = State(0.0), State(0.0), State(0.0)
    for target, key, name in [(lox_target, "ox_mass_flow", "OTV Flow Target"),
                              (fuel_target, "fuel_mass_flow", "FTV Flow Target"),
                              (thrust_target, "thrust", "Thrust Target")]:
        times, values = prop.target_schedule(targets, key, config["thrust_profile"]["interpolation"])
        Sequence(name, net, target=target, times=times, values=values)

    # ---- Reservoirs and main feeds ---- #
    PC, PT = State(D.AMBIENT_PRESSURE), State(D.AMBIENT_PRESSURE)   # main and torch chamber pressure [Pa]
    tanks, feeds = {}, {}
    for label, side, tank, elements, target in [("LOX", "oxidizer", D.LOX_TANK, hardware["oxidizer_feed"], lox_target),
                                                 ("Fuel", "fuel", D.FUEL_TANK, hardware["fuel_feed"], fuel_target)]:
        spec = reservoirs[side]
        liquid = FunctionLookup(f"{spec['name']} Liquid", net, prop.liquid_properties, propellant=tank.propellant,
                                pressure=spec["pressure_psia"] * psi, temperature=spec["temperature_k"],
                                outputs=prop.LIQUID_PROPERTIES)
        tanks[label] = {"liquid": liquid, "liquid_node": SimpleNamespace(pressure=State(spec["pressure_psia"] * psi))}
        feeds[label] = prop.build_liquid_feed(net, label, elements, tanks[label], PC, target)

    def node_pressure(label, name):
        return dict(feeds[label]["nodes"])[name]

    # ---- Torch legs, chamber and throat ---- #
    torch_flows = {}
    for label, key in (("LOX", "torch_oxidizer_leg"), ("Fuel", "torch_fuel_leg")):
        leg = hardware[key]
        source = node_pressure(label, leg["source"])
        name = "Torch Ox" if label == "LOX" else "Torch Fuel"
        torch_flows[label] = build_torch_leg(net, name, leg["elements"], source, tanks[label]["liquid"], PT,
                                             prop.track_element)
        junction, next_flow = feeds[label]["junctions"][leg["source"]]
        junction.mass_flow_out = next_flow + torch_flows[label]
    floor = D.MIXTURE_RATIO_FLOW_FLOOR
    torch_mr = ((torch_flows["LOX"] + D.MIXTURE_RATIO * floor) / (torch_flows["Fuel"] + floor)).clip(
        D.CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT, D.CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT)
    torch_gas = Map.from_hdf5("Torch Gas", net, D.CHAMBER_MAP_FILENAME, group="chamber",
                              inputs={"chamber_pressure": PT, "mixture_ratio": torch_mr}, extrapolate=True)

    # ---- Main chamber ---- #
    ox_in = feeds["LOX"]["chamber_inflow"] + torch_flows["LOX"]
    fuel_in = feeds["Fuel"]["chamber_inflow"] + torch_flows["Fuel"]
    main_mr = ((ox_in + D.MIXTURE_RATIO * floor) / (fuel_in + floor)).clip(
        D.CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT, D.CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT)
    engine = EnginePerformance("Main Chamber", net, chamber_pressure=PC, mixture_ratio=main_mr,
                               sampled_efficiency=True, chamber_map_filename=D.CHAMBER_MAP_FILENAME)
    throat = FunctionLookup("Torch Throat", net, torch_throat_flow, p_up=PT, p_down=PC,
                            gamma_up=torch_gas.gamma, r_up=torch_gas.gas_constant, t_up=torch_gas.temperature,
                            gamma_down=engine["gas"].gamma, r_down=engine["gas"].gas_constant,
                            t_down=engine["gas"].temperature, area=torch["throat_area_in2"] * 0.0254**2,
                            band=torch["near_zero_dp_psi"] * psi, outputs=("mass_flow",))
    torch_out = throat.mass_flow / torch["cstar_efficiency"]
    Volume(torch["chamber_name"], net, pressure=PT, mass_flow_in=torch_flows["LOX"] + torch_flows["Fuel"],
           mass_flow_out=torch_out)
    # The torch products arrive through the throat; at steady state that equals the torch
    # inflow, so the main chamber counts the torch ox and fuel directly in its flow and MR.
    chamber = Volume("Main Chamber", net, pressure=PC, mass_flow_in=feeds["LOX"]["chamber_inflow"]
                     + feeds["Fuel"]["chamber_inflow"] + torch_out)
    chamber.mass_flow_out = engine["mass_flow"]

    # ---- Purge system ---- #
    supply = config["purge_supply"]
    gas = Lookup(f"{supply['name']} Gas", net, Fluid, D.PRESSURANT, pressure=supply["pressure_psia"] * psi,
                 temperature=supply["temperature_k"])
    nodes = {supply["name"]: {"pressure": State(supply["pressure_psia"] * psi), "fluid": gas, "gas": True,
                              "store": supply["name"]}}
    sinks = {"ambient": D.AMBIENT_PRESSURE, torch["chamber_name"]: PT, "Main Chamber": PC}
    for label in feeds:
        sinks.update(dict(feeds[label]["nodes"]))
    prop.build_branches(net, [hardware["regulator"]], hardware["branches"], nodes, sinks)

    # ---- Tracks ---- #
    lbf = pound_force
    net.track("Chamber Pressure [psia]", PC / psi)
    net.track("Torch Chamber Pressure [psia]", PT / psi)
    net.track("Thrust [lbf]", engine["thrust"] / lbf)
    net.track("Target Thrust [lbf]", thrust_target / lbf)
    net.track("Main LOX Mass Flow [kg/s]", feeds["LOX"]["chamber_inflow"])
    net.track("Main Fuel Mass Flow [kg/s]", feeds["Fuel"]["chamber_inflow"])
    net.track("Torch LOX Mass Flow [kg/s]", torch_flows["LOX"])
    net.track("Torch Fuel Mass Flow [kg/s]", torch_flows["Fuel"])
    net.track("LOX Into Chamber [kg/s]", ox_in)
    net.track("Fuel Into Chamber [kg/s]", fuel_in)
    net.track("LOX Reservoir Outflow [kg/s]", feeds["LOX"]["tank_outflow"])
    net.track("Fuel Reservoir Outflow [kg/s]", feeds["Fuel"]["tank_outflow"])
    net.track("Torch Throat Mass Flow [kg/s]", torch_out)
    net.track("Mixture Ratio", main_mr)
    net.track("Torch Mixture Ratio", torch_mr)
    net.track("Torch Temperature [K]", torch_gas.temperature)
    net.track("Chamber Temperature [K]", engine["gas"].temperature)
    net.track("Nozzle Mass Flow [kg/s]", engine["mass_flow"])
    net.track("C* Efficiency", engine["cstar_efficiency"])
    net.track("Cf Efficiency", engine["cf_efficiency"])

    # ---- Solve ---- #
    dt = config["simulation"]["maximum_step_s"]
    filename = output / "Hopper_Hotfire.h5"
    print(f"Network ready in {time.perf_counter() - start:.1f} s; solving {t_final:g} s at <= {dt:g} s steps.",
          flush=True)
    transient_solver(net, max_evaluations=D.PROP_MAX_NFEV).solve(
        dt=dt, t_final=t_final, save_dt=config["simulation"]["save_interval_s"], filename=str(filename),
        rtol=D.PROP_RTOL, jacobian_method=D.PROP_JACOBIAN_METHOD, verbose=False, statistics=True)
    print(f"Solved in {time.perf_counter() - start:.1f} s.", flush=True)
    return filename, targets, timeline


# ---- Results ---- #
COLUMNS = {   # CSV column -> track
    "Thrust [lbf]": "Thrust [lbf]", "Target Thrust [lbf]": "Target Thrust [lbf]",
    "Chamber Pressure [psia]": "Chamber Pressure [psia]", "Torch Chamber Pressure [psia]": "Torch Chamber Pressure [psia]",
    "LOX Mass Flow, main injector [kg/s]": "Main LOX Mass Flow [kg/s]",
    "Fuel Mass Flow, main injector [kg/s]": "Main Fuel Mass Flow [kg/s]",
    "LOX Mass Flow, torch [kg/s]": "Torch LOX Mass Flow [kg/s]",
    "Fuel Mass Flow, torch [kg/s]": "Torch Fuel Mass Flow [kg/s]",
    "LOX Target [kg/s]": "OTV Target Mass Flow [kg/s]", "Fuel Target [kg/s]": "FTV Target Mass Flow [kg/s]",
    "OTV Open [%]": "OTV Opening", "FTV Open [%]": "FTV Opening",
    "LOX Injector dP [psi]": "LOX Injector dP [psi]", "Fuel Injector dP [psi]": "Fuel Injector dP [psi]",
    "Mixture Ratio": "Mixture Ratio", "Torch Mixture Ratio": "Torch Mixture Ratio",
    "LOX Purge [kg/s]": "LOX Purge Mass Flow [kg/s]", "Fuel Purge [kg/s]": "Fuel Purge Mass Flow [kg/s]",
    "Torch Purge [kg/s]": "Torch Purge Mass Flow [kg/s]", "Purge Relief [kg/s]": "Purge Relief Mass Flow [kg/s]",
    "Purge Reg Outlet [psia]": "Purge Reg Outlet Pressure [psia]",
}


def read_history(filename):
    import h5py
    import numpy as np
    from Backend.Hopper_Plots import _safe_name
    with h5py.File(filename, "r") as h5:
        network = next(name for name, group in h5.items() if isinstance(group, h5py.Group) and "transient" in group)
        base = h5[f"{network}/transient/runs/base"]
        history = {"Time [s]": np.asarray(base["time"], dtype=float)}
        for column, track in COLUMNS.items():
            history[column] = np.asarray(base["tracks"][_safe_name(track)], dtype=float)
    for valve in ("OTV", "FTV"):
        history[f"{valve} Open [%]"] = 100 * history[f"{valve} Open [%]"]
    return history


def summarize(history, timeline, config):
    import numpy as np
    t = history["Time [s]"]

    def at(column, when):
        return float(np.interp(when, t, history[column]))

    end = timeline["profile_end_s"]
    burning = (t >= timeline["profile_start_s"]) & (t <= end)
    pc = history["Chamber Pressure [psia]"]
    stiffness = {k: float(np.min(history[f"{k} Injector dP [psi]"][burning] / pc[burning]))
                 for k in ("LOX", "Fuel")} if burning.any() else {}
    dt = np.diff(t)

    def total(column):
        values = history[column]
        return float(np.sum(0.5 * (values[1:] + values[:-1]) * dt))

    return {
        "timeline_s": timeline,
        "at_end_of_ramp": {column: at(column, end) for column in (
            "Thrust [lbf]", "Target Thrust [lbf]", "Chamber Pressure [psia]", "LOX Mass Flow, main injector [kg/s]",
            "Fuel Mass Flow, main injector [kg/s]", "OTV Open [%]", "FTV Open [%]", "LOX Injector dP [psi]",
            "Fuel Injector dP [psi]", "Mixture Ratio")},
        "minimum_injector_stiffness_during_profile": stiffness,
        "peak": {column: float(np.max(history[column])) for column in (
            "Thrust [lbf]", "Chamber Pressure [psia]", "Torch Chamber Pressure [psia]", "OTV Open [%]", "FTV Open [%]")},
        "torch_mid_window": {column: at(column, 0.5 * (timeline["torch_open_s"] + min(timeline["main_fuel_opens_s"], timeline["main_ox_opens_s"])))
                             for column in ("Torch Chamber Pressure [psia]", "LOX Mass Flow, torch [kg/s]",
                                            "Fuel Mass Flow, torch [kg/s]", "Torch Mixture Ratio")},
        "propellant_used_kg": {"LOX": total("LOX Mass Flow, main injector [kg/s]") + total("LOX Mass Flow, torch [kg/s]"),
                               "Fuel": total("Fuel Mass Flow, main injector [kg/s]") + total("Fuel Mass Flow, torch [kg/s]")},
        "purge_gas_used_kg": {name: total(f"{name} [kg/s]") for name in
                              ("LOX Purge", "Fuel Purge", "Torch Purge", "Purge Relief")},
        "reservoir_pressures_psia": {k: v["pressure_psia"] for k, v in config["reservoirs"].items()},
        "assumptions": ASSUMPTIONS,
    }


PROFILES = {   # file -> graphed columns (the plots in plot())
    "mass_flow": ["LOX Mass Flow, main injector [kg/s]", "Fuel Mass Flow, main injector [kg/s]",
                  "LOX Mass Flow, torch [kg/s]", "Fuel Mass Flow, torch [kg/s]"],
    "throttle_valves": ["OTV Open [%]", "FTV Open [%]"],
    "thrust": ["Thrust [lbf]"],
    "chamber_pressure": ["Chamber Pressure [psia]"],
    "injector_dp": ["LOX Injector dP [psi]", "Fuel Injector dP [psi]"],
}
LBM_PER_KG = 2.2046226218


def export_profiles(history, output):
    """One CSV per graph with its full time history (mass flows also in lbm/s),
    and hotfire_maxima.csv: the peak of every graphed quantity and when it occurs."""
    import numpy as np
    folder = output / "profiles"
    folder.mkdir(exist_ok=True)
    t = history["Time [s]"]
    maxima = []
    for key, columns in PROFILES.items():
        table = {"Time [s]": t}
        for column in columns:
            table[column] = history[column]
            if column.endswith("[kg/s]"):
                table[column.replace("[kg/s]", "[lbm/s]")] = history[column] * LBM_PER_KG
        with (folder / f"{key}.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(list(table))
            writer.writerows(zip(*table.values()))
        for column in table:
            if column == "Time [s]":
                continue
            i = int(np.argmax(table[column]))
            maxima.append({"graph": key, "quantity": column, "maximum": float(table[column][i]),
                           "time [s]": float(t[i])})
    with (output / "hotfire_maxima.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(maxima[0]))
        writer.writeheader()
        writer.writerows(maxima)
    return maxima


def plot(history, timeline, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = history["Time [s]"]
    lb = 2.2046226218
    panels = [
        ("mass_flow", "Mass flow [kg/s]", [
            ("LOX Mass Flow, main injector [kg/s]", "LOX (main injector)", "C0", "-"),
            ("LOX Target [kg/s]", "LOX target", "C0", ":"),
            ("Fuel Mass Flow, main injector [kg/s]", "Fuel (main injector)", "C1", "-"),
            ("Fuel Target [kg/s]", "Fuel target", "C1", ":"),
            ("LOX Mass Flow, torch [kg/s]", "LOX (torch)", "C0", "--"),
            ("Fuel Mass Flow, torch [kg/s]", "Fuel (torch)", "C1", "--")]),
        ("throttle_valves", "Throttle valve open [%]", [
            ("OTV Open [%]", "OTV", "C0", "-"), ("FTV Open [%]", "FTV", "C1", "-")]),
        ("thrust", "Thrust [lbf]", [
            ("Thrust [lbf]", "Thrust", "k", "-"), ("Target Thrust [lbf]", "Target", "C3", ":")]),
        ("chamber_pressure", "Chamber pressure [psia]", [
            ("Chamber Pressure [psia]", "Main chamber Pc", "k", "-")]),
        ("injector_dp", "Injector dP [psi]", [
            ("LOX Injector dP [psi]", "LOX (OINJ)", "C0", "-"), ("Fuel Injector dP [psi]", "Fuel (FINJ)", "C1", "-")]),
    ]

    def decorate(ax):
        ax.axvspan(timeline["torch_open_s"], timeline["torch_close_s"], color="C4", alpha=0.07, lw=0)
        ax.axvspan(timeline["profile_start_s"], timeline["profile_end_s"], color="0.5", alpha=0.07, lw=0)
        ax.axvspan(timeline["main_purges_open_s"], timeline["main_purges_close_s"], color="C2", alpha=0.08, lw=0)
        ax.grid(alpha=0.3)

    def draw(ax, key, label, lines):
        for column, name, color, style in lines:
            ax.plot(t, history[column], color=color, ls=style, lw=1.6 if style == "-" else 1.2, label=name)
        ax.set_ylabel(label)
        ax.legend(fontsize=8, loc="best")
        decorate(ax)
        if key == "mass_flow":
            second = ax.secondary_yaxis("right", functions=(lambda x: x * lb, lambda x: x / lb))
            second.set_ylabel("[lbm/s]")

    figure, axes = plt.subplots(len(panels), 1, figsize=(10, 15), sharex=True)
    for ax, (key, label, lines) in zip(axes, panels):
        draw(ax, key, label, lines)
    axes[0].set_title("Test-stand hotfire  (purple: torch on, grey: thrust profile, green: OPV/FPV purge)")
    axes[-1].set_xlabel("Time [s]")
    figure.tight_layout()
    figure.savefig(output / "hotfire_overview.png", dpi=140)
    plt.close(figure)
    for key, label, lines in panels:
        figure, ax = plt.subplots(figsize=(9, 4.5))
        draw(ax, key, label, lines)
        ax.set_xlabel("Time [s]")
        figure.tight_layout()
        figure.savefig(output / f"hotfire_{key}.png", dpi=140)
        plt.close(figure)


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

    if args.check_inputs:
        run(config, output, check_inputs=True)
        return
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(config_path, output / "hotfire_config.json")
    import HOPPER_DATA as D
    shutil.copy2(D.CONFIG_PATH, output / "hopper_config.json")
    filename, targets, timeline = run(config, output)
    history = read_history(filename)
    with (output / "hotfire_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(list(history))
        writer.writerows(zip(*history.values()))
    with (output / "thrust_targets.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["profile time [s]", "thrust [lbf]", "Pc [psia]", "LOX [kg/s]", "fuel [kg/s]"])
        for p in targets:
            writer.writerow([p["time"], p["thrust"] / 4.4482216152605, p["chamber_pressure"] / 6894.757293168361,
                             p["ox_mass_flow"], p["fuel_mass_flow"]])
    summary = summarize(history, timeline, config)
    summary["elapsed_s"] = time.perf_counter() - started
    (output / "hotfire_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    plot(history, timeline, output)
    export_profiles(history, output)
    print(json.dumps({k: summary[k] for k in ("at_end_of_ramp", "minimum_injector_stiffness_during_profile",
                                              "torch_mid_window", "peak")}, indent=2))
    print(f"Results in {output}")


if __name__ == "__main__":
    main()

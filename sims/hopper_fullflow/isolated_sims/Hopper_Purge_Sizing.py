"""Purge sizing: the smallest CdA of the sized purge elements (the check valves)
that still meets the purge requirement, and the regulator set it then needs.

Inputs are the purge section of hopper_config.json:
    sequence                    branches purged, each opening at its delay after the shutdown
                                request for its duration; overlapping purges (and the torch
                                purge, open throughout) share the regulator's capacity
    required_fraction_of_peak_flow, back_pressure_psia
                                each purge's gas volume flow at back pressure must be at
                                least this fraction of its propellant's peak volume flow,
                                for the whole time it is open
    sizing.elements             elements given one shared CdA, the unknown (CV-4, CV-5)
    sizing.regulator_min_dp_psi margin the regulator needs above its set to regulate
    sizing.max_regulator_set_psig   optional ceiling on the set (e.g. below ORV-1)

From the last full run it takes the COPV state at the shutdown request,
each propellant's peak volume flow, and the feed-line pressure at that peak (the back
pressure, unless back_pressure_psia sets one).

Model
-----
COPV: the gas left in the bottles expands isentropically as it is drawn off (an
adiabatic, well-mixed blowdown), stepped every sizing.time_step_s. Regulator: outlet =
ambient + set (+ its supply-pressure effect) while the COPV is at least that + min dP;
its CdA caps the combined flow of every open purge (below set when exceeded). Purge line: its restrictions in series (FullFlow
SeriesCdA rule), compressible orifice flow with the real-gas R and gamma, gas throttled
at constant enthalpy (ThermoProp), as in Hopper_Purge.py.

For a trial CdA, the required set is the lowest that keeps every purge at or above the
requirement. A set that is too high draws the COPV below set + min dP before the
sequence ends (the regulator drops out), so a CdA is feasible only if the lowest set
that meets the flow requirement still holds regulation to the end. The smallest CdA is
found by bisection, the set by bisection inside it.

Outputs in PROP_RESULTS_DIRECTORY:
    purge_sizing.csv            required regulator set vs CdA of the sized elements
    purge_sizing_sequence.csv   the purge sequence at the minimum CdA, step by step
    purge_sizing.png            both, plotted

    python isolated_sims/Hopper_Purge_Sizing.py [results.h5] [options]
    options (override the config for this run only):
      --durations 2 2            seconds each purge stays open, in config order
      --delays 0 0.1             seconds after the shutdown request each purge opens
      --fraction 1.0             required fraction of peak propellant volume flow
      --back-pressure 400        back pressure [psia] for every purge (default: from the run)
      --elements CV-4 CV-5       elements whose shared CdA is solved for
      --min-dp 0                 regulator margin above set [psi]
      --max-set 800              ceiling on the regulator set [psig]
"""

import argparse
import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

if not __package__:   # run as a script: make HOPPER_DATA and the Backend package importable
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import HOPPER_DATA as D
from scipy.constants import liter, milli, psi
from thermoprop import Fluid

from Backend.Hopper_Elements import ReliefValve
from Backend.Hopper_Export import LoadTracks, _track
from isolated_sims.Hopper_Purge import orifice_mass_flow, purge_start_time

FEED_TRACKS = {"LOX": ("Ox Mass Flow [kg/s]", "LOX Liquid Density [kg/m3]"),
               "Fuel": ("Fuel Mass Flow [kg/s]", "Fuel Liquid Density [kg/m3]")}


# ---- Inputs from the config and the run ---- #
def _branch(name):
    return next(branch for branch in D.BRANCHES if branch.name == name)


def _regulator(branch):
    """The regulator whose outlet feeds this branch, or None."""
    return next((r for r in D.REGULATORS if r.outlet_node == branch.source), None)


def _side(branch):
    """"LOX" or "Fuel": the feed line the purge discharges into."""
    for label, feed in (("LOX", D.LOX_FEED), ("Fuel", D.FUEL_FEED)):
        if branch.sink in [element.outlet_node for element in feed]:
            return label
    raise ValueError(f"{branch.name}: sink {branch.sink!r} is not a LOX or fuel feed junction, so it has no "
                     "propellant flow to match.")


def run_state(filename=None):
    """COPV state when the first purge opens, and each purged line's peak propellant
    volume flow [m^3/s] and the feed pressure [Pa] at that moment."""
    filename = filename or D.PROP_FILENAME
    time, tracks = LoadTracks(filename)
    start = purge_start_time(filename, time)
    total = _track(tracks, "Total Mass Flow [kg/s]")
    burning = (total > 0.05 * total.max()) & (time <= start)
    legs = {}
    for name, _, _ in D.PURGE_SEQUENCE:
        branch = _branch(name)
        flow_track, density_track = FEED_TRACKS[_side(branch)]
        volume_flow = _track(tracks, flow_track) / _track(tracks, density_track)
        peak = int(np.argmax(np.where(burning, volume_flow, -1.0)))
        legs[name] = {"peak_flow": float(volume_flow[peak]),
                      "feed_pressure": float(_track(tracks, f"{branch.sink} Pressure [Pa]")[peak])}
    return {"copv_pressure": float(np.interp(start, time, _track(tracks, "COPV Pressure [psia]"))) * psi,
            "copv_temperature": float(np.interp(start, time, _track(tracks, "COPV Temperature [K]"))),
            "purge_start": start, "legs": legs}


def path_cda(branch, sized_cda):
    """Series CdA [m^2] of a branch with the sized elements at sized_cda (None: as configured)."""
    return sum(1.0 / (sized_cda if sized_cda is not None and element.name in D.PURGE_SIZING_ELEMENTS
                      else element.cda_value) ** 2
               for element in branch.elements) ** -0.5


# ---- One purge sequence ---- #
def _torch_purge_cda():
    """Series CdA [m^2] of the torch purge (TFPV path and the shared torch fuel injector), or None."""
    if not D.TORCH_LEGS:
        return None
    purge = D.TORCH_LEGS["purge"]
    return sum(1.0 / e.cda_value ** 2 for e in purge["elements"] + purge["shared"]) ** -0.5


def _solve_outlet(target, back_floor, demand, capacity):
    """Regulator outlet [Pa]: the target if the regulator can pass the demand there,
    else the pressure at which its capacity equals the demand (bisection)."""
    if capacity is None or demand(target) <= capacity(target):
        return target
    low, high = back_floor, target
    for _ in range(60):
        middle = 0.5 * (low + high)
        if demand(middle) > capacity(middle):
            high = middle
        else:
            low = middle
    return low


def simulate(state, set_gauge, sized_cda, record=False, regulator_cda=None, required_mass_flow=None):
    """Blow the COPV down through the purge sequence at one regulator set [Pa gauge].

    Purges open at their open_after_shutdown_s and may overlap; every open purge,
    and the torch purge (open throughout, discharging to ambient: the most it can
    draw), shares the regulator. The regulator holds set (plus its supply-pressure
    effect) while its capacity covers the combined demand; past that its outlet
    falls to where capacity = demand.

    regulator_cda [m^2] replaces the regulator's configured CdA (regulator sizing).
    required_mass_flow {purge: kg/s} replaces the volume-flow requirement with a mass
    flow each purge must reach.

    Returns {"status": "ok" | "short" | "dropout", "min_ratio", "copv_end_pressure",
    "gas_used", "rows"}. "short": some purge fell below the requirement. "dropout":
    the COPV fell below the target + min dP (the set is too high for the gas left),
    checked before "short" because a higher set cannot fix it."""
    volume = D.COPV_VOLUME
    copv = Fluid(D.PRESSURANT, pressure=state["copv_pressure"], temperature=state["copv_temperature"])
    entropy, mass = copv.entropy, copv.density * volume
    start_mass, min_ratio, dropped, rows = mass, math.inf, False, []
    windows = [(name, opened, opened + duration) for name, opened, duration in D.PURGE_SEQUENCE]
    regulator = _regulator(_branch(windows[0][0]))
    if any(_regulator(_branch(name)) is not regulator for name, _, _ in windows):
        raise ValueError("purge sizing: every purge must draw from the same regulator.")
    reference = (regulator.reference_inlet_pressure or D.PRESSURANT_INITIAL_PRESSURE) if regulator else 0.0
    torch_cda = _torch_purge_cda()
    legs = {}
    for name, _, _ in windows:
        leg = state["legs"][name]
        legs[name] = (path_cda(_branch(name), sized_cda),
                      D.PURGE_BACK_PRESSURE if D.PURGE_BACK_PRESSURE is not None else leg["feed_pressure"],
                      D.PURGE_REQUIRED_FRACTION * leg["peak_flow"])
    end = max(closed for _, _, closed in windows)
    steps = max(1, math.ceil(end / D.PURGE_SIZING_DT))
    dt = end / steps
    for k in range(steps + 1):   # the last point is the state as the last valve closes
        t = k * dt
        copv = Fluid(D.PRESSURANT, density=mass / volume, entropy=entropy)
        active = [name for name, opened, closed in windows if opened <= t + 1e-9 and t <= closed + 1e-9]
        if regulator:
            target = D.AMBIENT_PRESSURE + set_gauge + regulator.supply_pressure_effect * (reference - copv.pressure)
            held = min(target, copv.pressure - D.PURGE_REGULATOR_MIN_DP)
            dropped = dropped or (active and held < target - 1e-3)
        else:
            held = copv.pressure

        def flows(outlet):
            gas = Fluid(D.PRESSURANT, pressure=outlet, enthalpy=copv.enthalpy)   # throttled at constant h
            out = {name: orifice_mass_flow(legs[name][0], outlet, gas.temperature, legs[name][1], gas.gamma,
                                           gas.gas_constant) for name in active}
            if torch_cda:
                out["Torch Purge"] = orifice_mass_flow(torch_cda, outlet, gas.temperature, D.AMBIENT_PRESSURE,
                                                       gas.gamma, gas.gas_constant)
            return out

        capacity, reg_cda = None, regulator_cda or (regulator.cda_value if regulator is not None else None)
        if reg_cda:
            capacity = lambda outlet: orifice_mass_flow(reg_cda, copv.pressure, copv.temperature, outlet,
                                                        copv.gamma, copv.gas_constant)
        floor = max([legs[name][1] for name in active] + [D.AMBIENT_PRESSURE])
        outlet = _solve_outlet(held, floor, lambda pressure: sum(flows(pressure).values()), capacity) \
            if (active or torch_cda) else held
        flow = flows(outlet)
        for name in active:
            back = legs[name][1]
            volume_flow = flow[name] / Fluid(D.PRESSURANT, pressure=back, enthalpy=copv.enthalpy).density
            ratio = (flow[name] / required_mass_flow[name] if required_mass_flow
                     else volume_flow / legs[name][2])
            min_ratio = min(min_ratio, ratio)
            if record:
                rows.append({"Time [s]": round(state["purge_start"] + t, 3), "Purge": name,
                             "COPV [psia]": round(copv.pressure / psi, 1), "COPV [K]": round(copv.temperature, 1),
                             "Regulator outlet [psia]": round(outlet / psi, 1),
                             "Back pressure [psia]": round(back / psi, 1),
                             "Purge mass flow [g/s]": round(1000 * flow[name], 2),
                             "Torch purge [g/s]": round(1000 * flow.get("Torch Purge", 0.0), 2),
                             "Purge volume flow [L/s]": round(volume_flow / liter, 4),
                             "Required [L/s]": round(legs[name][2] / liter, 4), "Purge / required": round(ratio, 4)})
        if k < steps:
            mass -= sum(flow.values()) * dt
    status = "dropout" if dropped else ("short" if min_ratio < 1.0 else "ok")
    return {"status": status, "min_ratio": min_ratio, "copv_end_pressure": copv.pressure,
            "gas_used": start_mass - mass, "rows": rows}


def required_set(state, sized_cda, tolerance=0.2 * psi):
    """Lowest regulator set [Pa gauge] meeting the requirement at this CdA, and whether
    the regulator holds it to the end. (None, result) if no set can meet it."""
    backs = [D.PURGE_BACK_PRESSURE if D.PURGE_BACK_PRESSURE is not None else leg["feed_pressure"]
             for leg in state["legs"].values()]
    low = max(backs) - D.AMBIENT_PRESSURE
    high = state["copv_pressure"] - D.PURGE_REGULATOR_MIN_DP - D.AMBIENT_PRESSURE
    if high <= low or simulate(state, high, sized_cda)["status"] == "short":
        return None, simulate(state, max(high, low), sized_cda)
    while high - low > tolerance:   # short below the answer, not short above it
        middle = 0.5 * (low + high)
        if simulate(state, middle, sized_cda)["status"] == "short":
            low = middle
        else:
            high = middle
    return high, simulate(state, high, sized_cda, record=True)


def feasible(state, sized_cda):
    set_gauge, result = required_set(state, sized_cda)
    ok = set_gauge is not None and result["status"] == "ok"
    if ok and D.PURGE_MAX_REGULATOR_SET is not None:
        ok = set_gauge <= D.PURGE_MAX_REGULATOR_SET
    return ok, set_gauge, result


def minimum_cda(state, relative_tolerance=0.002):
    """(smallest feasible CdA [m^2], its set [Pa gauge], its sequence), or Nones if even the
    top of the range fails. If the bottom of the range already works, that is returned
    with the flag at_range_limit."""
    low, high = D.PURGE_SIZING_CDA_RANGE
    ok, set_gauge, result = feasible(state, high)
    if not ok:
        return {"cda": None, "set": None, "result": result, "at_range_limit": False}
    ok_low, set_low, result_low = feasible(state, low)
    if ok_low:
        return {"cda": low, "set": set_low, "result": result_low, "at_range_limit": True}
    best = (high, set_gauge, result)
    while high / low - 1.0 > relative_tolerance:
        middle = math.sqrt(low * high)
        ok, set_gauge, result = feasible(state, middle)
        if ok:
            high, best = middle, (middle, set_gauge, result)
        else:
            low = middle
    return {"cda": best[0], "set": best[1], "result": best[2], "at_range_limit": False}


# ---- Report ---- #
def _relief_crack():
    """Crack pressure [Pa gauge] of a relief on the regulator outlet, or None."""
    outlets = {r.outlet_node for r in D.REGULATORS}
    for branch in D.BRANCHES:
        if branch.source in outlets:
            for element in branch.elements:
                if isinstance(element, ReliefValve):
                    return element.name, element.crack_pressure
    return None, None


def purge_sizing(filename=None, directory=None):
    """Solve, write purge_sizing.csv / _sequence.csv / .png; returns {metric: value}."""
    import matplotlib.pyplot as plt

    directory = directory or D.PROP_RESULTS_DIRECTORY
    os.makedirs(directory, exist_ok=True)
    state = run_state(filename)
    answer = minimum_cda(state)
    sized = [e for name, _, _ in D.PURGE_SEQUENCE for e in _branch(name).elements if e.name in D.PURGE_SIZING_ELEMENTS]
    current = sized[0].cda_value
    relief, crack = _relief_crack()
    label = " / ".join(D.PURGE_SIZING_ELEMENTS)
    sequence = ", ".join(f"{name} {duration:g} s from T+{delay:g}" for name, delay, duration in D.PURGE_SEQUENCE)

    # Required set over a range of CdA, for the plot and the table.
    low, high = D.PURGE_SIZING_CDA_RANGE
    start = answer["cda"] or low
    sweep = sorted(set(np.geomspace(start, max(4 * start, 2 * current, start * 1.01), 6).tolist() + [current]))
    rows = []
    for cda in sweep:
        set_gauge, result = required_set(state, cda)
        rows.append({f"{label} CdA [mm2]": round(cda / milli**2, 3),
                     "Path CdA [mm2]": round(path_cda(_branch(D.PURGE_SEQUENCE[0][0]), cda) / milli**2, 3),
                     "Required regulator set [psig]": "none" if set_gauge is None else round(set_gauge / psi, 1),
                     "Regulator holds to the end": "n/a" if set_gauge is None else ("yes" if result["status"] == "ok" else "NO"),
                     "COPV after purges [psia]": round(result["copv_end_pressure"] / psi, 1),
                     "Gas used [g]": round(1000 * result["gas_used"], 1)})
    with open(os.path.join(directory, "purge_sizing.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    design = answer["result"]["rows"] if answer["cda"] else []
    if design:
        with open(os.path.join(directory, "purge_sizing_sequence.csv"), "w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=list(design[0]))
            writer.writeheader()
            writer.writerows(design)

    summary = {"Purge sequence": sequence,
               "Purge requirement": f">= {100 * D.PURGE_REQUIRED_FRACTION:g} % of peak propellant volume flow",
               "COPV when purges open [psia]": round(state["copv_pressure"] / psi, 1),
               "COPV temperature when purges open [K]": round(state["copv_temperature"], 1)}
    for name, leg in state["legs"].items():
        back = D.PURGE_BACK_PRESSURE if D.PURGE_BACK_PRESSURE is not None else leg["feed_pressure"]
        summary[f"{name} required flow [L/s] at [psia]"] = (f"{D.PURGE_REQUIRED_FRACTION * leg['peak_flow'] / liter:.3f} "
                                                            f"at {back / psi:.1f}")
    if answer["cda"] is None:
        summary[f"Minimum {label} CdA [mm2]"] = (f"none up to {high / milli**2:g} mm2: even at the highest set the "
                                                 "COPV allows, a purge falls short. Usually the regulator's capacity "
                                                 "(see regulator_sizing.csv / Hopper_Regulator_Sizing.py)")
    else:
        bound = " (bottom of the search range)" if answer["at_range_limit"] else ""
        summary[f"Minimum {label} CdA [mm2]"] = f"{answer['cda'] / milli**2:.3f}{bound}"
        summary["Purge path CdA at minimum [mm2]"] = round(path_cda(_branch(D.PURGE_SEQUENCE[0][0]), answer["cda"]) / milli**2, 3)
        summary["Regulator set at minimum CdA [psig]"] = round(answer["set"] / psi, 1)
        summary["COPV after purges at minimum CdA [psia]"] = round(answer["result"]["copv_end_pressure"] / psi, 1)
        summary["Purge gas used at minimum CdA [g]"] = round(1000 * answer["result"]["gas_used"], 1)
    current_set, current_result = required_set(state, current)
    summary[f"Regulator set needed at current {label} {current / milli**2:g} mm2 [psig]"] = (
        "none" if current_set is None else round(current_set / psi, 1))
    if crack is not None and answer["set"] is not None and answer["set"] >= crack:
        summary["Note"] = (f"the minimum-CdA set is at or above {relief}'s {crack / psi:.0f} psig crack; raise {relief} "
                           "or set purge.sizing.max_regulator_set_psig below it")

    # ---- Figure ---- #
    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    ax = axes[0]
    sets = [r["Required regulator set [psig]"] for r in rows]
    xs = [r[f"{label} CdA [mm2]"] for r in rows]
    ax.plot([x for x, s in zip(xs, sets) if s != "none"], [s for s in sets if s != "none"], "o-", label="required regulator set")
    ax.plot(xs, [r["COPV after purges [psia]"] - D.AMBIENT_PRESSURE / psi - D.PURGE_REGULATOR_MIN_DP / psi for r in rows],
            "--", color="0.4", label="COPV after purges - margin (highest set it can hold)")
    if crack is not None:
        ax.axhline(crack / psi, color="tab:red", lw=1, label=f"{relief} crack")
    if answer["cda"]:
        ax.axvline(answer["cda"] / milli**2, color="k", lw=1, ls=":", label=f"minimum CdA {answer['cda'] / milli**2:.2f} mm²")
    ax.axvline(current / milli**2, color="tab:green", lw=1, ls=":", label=f"current {current / milli**2:g} mm²")
    ax.set_xlabel(f"{label} CdA [mm²]"); ax.set_ylabel("Pressure [psig]"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax.set_title("(a) Regulator set each check-valve CdA needs")
    ax = axes[1]
    if design:
        t = [r["Time [s]"] for r in design]
        ax.plot(t, [r["Purge / required"] for r in design], label="purge / required volume flow")
        ax.axhline(1.0, color="k", lw=1)
        twin = ax.twinx()
        twin.plot(t, [r["COPV [psia]"] for r in design], color="0.6", label="COPV")
        twin.plot(t, [r["Regulator outlet [psia]"] for r in design], color="tab:orange", label="regulator outlet")
        twin.set_ylabel("Pressure [psia]"); twin.legend(fontsize=8, loc="center right")
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Purge / required"); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="lower left")
    ax.set_title(f"(b) At the minimum CdA: {sequence}")
    fig.suptitle("Purge sizing", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(directory, "purge_sizing.png"), dpi=120)
    plt.close(fig)
    return summary


def _apply_overrides(args):
    """Command-line values replace the config's for this run only."""
    sequence = list(D.PURGE_SEQUENCE)
    for option, index in (("durations", 2), ("delays", 1)):
        values = getattr(args, option)
        if values:
            if len(values) != len(sequence):
                raise SystemExit(f"--{option} needs {len(sequence)} values, one per purge in config order.")
            sequence = [tuple(v if i == index else leg[i] for i in range(3)) for leg, v in zip(sequence, values)]
    D.PURGE_SEQUENCE = sequence
    if args.fraction is not None:
        D.PURGE_REQUIRED_FRACTION = args.fraction
    if args.back_pressure is not None:
        D.PURGE_BACK_PRESSURE = args.back_pressure * psi
    if args.elements:
        D.PURGE_SIZING_ELEMENTS = args.elements
    if args.min_dp is not None:
        D.PURGE_REGULATOR_MIN_DP = args.min_dp * psi
    if args.max_set is not None:
        D.PURGE_MAX_REGULATOR_SET = args.max_set * psi


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    parser = argparse.ArgumentParser(description="Smallest purge check-valve CdA and the regulator set it needs.")
    parser.add_argument("results", nargs="?", help="results file (default: the latest run)")
    parser.add_argument("--durations", nargs="+", type=float)
    parser.add_argument("--delays", nargs="+", type=float)
    parser.add_argument("--fraction", type=float)
    parser.add_argument("--back-pressure", type=float)
    parser.add_argument("--elements", nargs="+")
    parser.add_argument("--min-dp", type=float)
    parser.add_argument("--max-set", type=float)
    arguments = parser.parse_args()
    _apply_overrides(arguments)
    for key, value in purge_sizing(arguments.results).items():
        print(f"{key:62s} {value}")
    print("wrote", os.path.join(D.PROP_RESULTS_DIRECTORY, "purge_sizing.png"))

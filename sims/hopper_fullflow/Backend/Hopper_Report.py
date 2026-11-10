"""Single-run reports: the numbers asked for after every Hopper run.

From a FullFlow results file and the thrust targets it flew:
    run_summary.csv   burn time, total impulse, propellant residuals (left in
                      the tanks at shutdown), COPV pressure / gas left, tank
                      pressure droop, injector stiffness, relief activity,
                      thrust-tracking and throttle-margin checks
    cda_table.csv     per thrust-profile point: target vs achieved thrust,
                      Pc, MR, flows, OTV/FTV CdA and % of their max CdA,
                      injector dP, tank pressures. This is the throttle-valve
                      CdA schedule, and the check that every thrust level is
                      reachable: a throttle at 100 % of max CdA has run out of
                      pressure and the thrust falls short.
    run_summary.txt / .json, propellant_residuals.csv
                      sustained demonstrated thrust, residuals, COPV temperature
                      limits and validity flags, from the metadata saved in the
                      results file
Vent times come from Hopper_Vent.py (vent_times.csv), LOX drain times from
Hopper_Drain.py (lox_drain_times.csv), and the tank pressure vs throttle CdA
sweep from Hopper_Tank_Pressure_Sweep.py (tank_pressure_sweep.csv / .png), and
the purge flow check from Hopper_Purge.py (purge_flow.csv / .png), and the purge
check-valve and regulator-set sizing from Hopper_Purge_Sizing.py (purge_sizing.csv / .png), and
the minimum regulator CdA / Cv from Hopper_Regulator_Sizing.py (regulator_sizing.csv / .png).

    python Backend/Hopper_Report.py [results.h5]
"""

import csv
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import psi

if not __package__:   # run as a script (python Backend/Hopper_Report.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Elements import ThrottleValve
from Backend.Hopper_Export import LoadTracks, _track
from Backend.Hopper_Inputs import LBF
from Backend.Hopper_Prop_System import burn_timeline, step_duration


# ---- Metadata saved with each run ---- #
def write_run_metadata(filename, targets, integration_dt):
    profile_start = D.PROFILE_START_TIME
    metadata = {
        "schema": 1, "profile": D.PROP_THRUST_CURVE_NAME,
        "profile_start": profile_start, "profile_end": profile_start + targets[-1]["time"],
        "integration_dt": integration_dt,
        "hold_time": D.MAX_THRUST_HOLD_TIME,
        "tracking_tolerance": D.MAX_THRUST_TRACKING_TOLERANCE,
        "map_pc_max_pa": D.CHAMBER_MAP_PRESSURE_UPPER_LIMIT,
        "map_mr_min": D.CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT,
        "map_mr_max": D.CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT,
        "copvs": [{"name": b.name, "minimum_temperature": b.min_temperature} for b in D.COPVS],
        "tanks": {"LOX": {"initial_mass": D.LOX_TANK.liquid_mass},
                  "Fuel": {"initial_mass": D.FUEL_TANK.liquid_mass}},
        "throttles": {label: [{"name": e.name, "max_cda": e.max_cda}
                                for e in feed if hasattr(e, "max_cda")]
                      for label, feed in [("LOX", D.LOX_FEED), ("Fuel", D.FUEL_FEED)]},
        "paths": {"ox": [e.name for e in D.LOX_FEED], "fuel": [e.name for e in D.FUEL_FEED],
                  "lox_pressurant": [e.name for e in D.LOX_PRESSURANT_LEG],
                  "fuel_pressurant": [e.name for e in D.FUEL_PRESSURANT_LEG]},
    }
    with h5py.File(filename, "a") as h5:
        h5.attrs["hopper_metadata"] = json.dumps(metadata)
        h5.attrs["hopper_config"] = json.dumps(getattr(D, 'RUN_CONFIG', D.CONFIG))


def read_metadata(filename):
    with h5py.File(filename, "r") as h5:
        raw = h5.attrs.get("hopper_metadata")
    return json.loads(raw) if raw else None


# ---- Sustained thrust, residuals, operating limits ---- #
def sustained_maximum(time, thrust, valid, duration):
    """Highest minimum thrust maintained across an entire duration.

    Conservative sampled definition: candidate intervals start/end at saved
    samples and span at least duration. Invalid samples break an interval.
    """
    time, thrust, valid = np.asarray(time), np.asarray(thrust), np.asarray(valid, dtype=bool)
    if duration <= 0 or len(time) < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("Sustained-thrust duration must be positive and saved times increasing.")
    best, interval = None, None
    for end in range(len(time)):
        start = int(np.searchsorted(time, time[end] - duration + 1e-10, side="right")) - 1
        if start < 0 or time[end] - time[start] < duration - 1e-9:
            continue
        if valid[start:end+1].all():
            value = float(np.min(thrust[start:end+1]))
            if best is None or value > best:
                best, interval = value, [float(time[start]), float(time[end])]
    return best, interval


def summarize(filename):
    time, tracks = LoadTracks(filename)
    metadata = read_metadata(filename)
    if metadata is None:
        raise ValueError("This history predates run metadata; rerun the simulation for a reliable maximum/residual report.")
    get = lambda name: _track(tracks, name)
    if metadata["hold_time"] <= 0:
        raise ValueError("MAX_THRUST_HOLD_TIME must be positive.")
    profile = (time >= metadata["profile_start"] - 1e-9) & (time <= metadata["profile_end"] + 1e-9)
    thrust, target = get("Thrust [N]"), get("Target Thrust [N]")
    valid = profile & np.isfinite(thrust)
    residuals = {}
    for label, tank in metadata["tanks"].items():
        mass = get(f"{label} Liquid Mass [kg]")
        valid &= mass > 0
        end_mass = float(np.interp(metadata["profile_end"], time, mass))
        residuals[label] = {"initial_kg": tank["initial_mass"], "profile_end_kg": end_mass,
                            "final_kg": float(mass[-1]), "consumed_kg": tank["initial_mass"] - float(mass[-1]),
                            "remaining_percent": 100 * float(mass[-1]) / tank["initial_mass"]}
    pc = get("Chamber Pressure [psia]") * 6894.757293168
    valid &= pc <= metadata["map_pc_max_pa"]
    # Use unclipped MR when available: the chamber map clips its input for numerical stability.
    fuel, ox = get("Fuel Mass Flow [kg/s]"), get("Ox Mass Flow [kg/s]")
    raw_mr = np.divide(ox, fuel, out=np.full_like(ox, np.inf), where=fuel > 1e-8)
    valid &= (raw_mr >= metadata["map_mr_min"]) & (raw_mr <= metadata["map_mr_max"])
    maximum, interval = sustained_maximum(time, thrust, valid, metadata["hold_time"])
    saturation = np.zeros(len(time), dtype=bool)
    for throttle_list in metadata["throttles"].values():
        for throttle in throttle_list:
            saturation |= get(f"{throttle['name']} CdA [mm2]") * 1e-6 >= 0.99 * throttle["max_cda"]
    shortfall = thrust < target * (1 - metadata["tracking_tolerance"])
    evidence, _ = sustained_maximum(time, np.ones_like(thrust), valid & saturation & shortfall,
                                    metadata["hold_time"])
    temperatures = []
    for bottle in metadata["copvs"]:
        values = get(f"{bottle['name']} Temperature [K]")
        limit = bottle["minimum_temperature"]
        crossing = np.flatnonzero(values < limit) if limit is not None else np.array([], dtype=int)
        temperatures.append({"name": bottle["name"], "minimum_observed_K": float(np.min(values)),
                             "limit_K": limit, "violated": bool(len(crossing)),
                             "first_crossing_s": float(time[crossing[0]]) if len(crossing) else None})
    return {"profile": metadata["profile"], "hold_time_s": metadata["hold_time"],
            "demonstrated_maximum_N": maximum, "demonstrated_maximum_lbf": maximum / LBF if maximum is not None else None,
            "interval_s": interval, "peak_thrust_N": float(np.max(thrust[profile])) if profile.any() else None,
            "capacity_shortfall_demonstrated": evidence is not None,
            "interpretation": ("Sustained valve saturation and thrust shortfall occurred; result applies to this run."
                               if evidence is not None else
                               "No sustained capacity shortfall demonstrated; achieved thrust is a lower bound, not the system maximum."),
            "residuals": residuals, "copv_temperatures": temperatures,
            "invalid_inventory": any(v["final_kg"] < 0 for v in residuals.values()),
            "profile_completed": bool(time[-1] >= metadata["profile_end"] - 1e-9),
            "excluded_profile_samples": int(np.count_nonzero(profile & ~valid))}


def export_summary(filename, output_dir):
    report = summarize(filename)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    (output / "run_summary.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    with (output / "propellant_residuals.csv").open("w", newline="") as stream:
        fields = ["propellant"] + list(next(iter(report["residuals"].values())))
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for name, values in report["residuals"].items():
            writer.writerow({"propellant": name, **values})
    thrust = report["demonstrated_maximum_N"]
    lines = [f"Profile: {report['profile']}"]
    if thrust is None:
        lines.append(f"No valid {report['hold_time_s']:g} s sustained-thrust interval was recorded.")
    else:
        lines.append(f"Max thrust demonstrated for >= {report['hold_time_s']:g} s: "
                     f"{thrust:.1f} N ({thrust / LBF:.1f} lbf)")
    lines.append(report["interpretation"])
    for label, values in report["residuals"].items():
        lines.append(f"{label} residual: {values['final_kg']:.4f} kg ({values['remaining_percent']:.2f} % initial); "
                     f"consumed {values['consumed_kg']:.4f} kg")
    for bottle in report["copv_temperatures"]:
        limit = "disabled" if bottle["limit_K"] is None else f"{bottle['limit_K']:g} K"
        lines.append(f"{bottle['name']}: minimum {bottle['minimum_observed_K']:.2f} K; limit {limit}"
                     + (f"; VIOLATED at {bottle['first_crossing_s']:.3f} s" if bottle["violated"] else ""))
    if report["invalid_inventory"]:
        lines.append("INVALID RUN: negative liquid inventory; results after depletion are outside the model.")
    if not report["profile_completed"]:
        lines.append("INCOMPLETE PROFILE: simulation stopped before the end of the commanded profile.")
    if report["excluded_profile_samples"]:
        lines.append(f"Excluded {report['excluded_profile_samples']} profile samples outside inventory/Pc/MR validity limits.")
    text = "\n".join(lines)
    print(text)
    (output / "run_summary.txt").write_text(text + "\n", encoding="utf-8")
    return report


# ---- Per-point CdA table and run summary ---- #
def _throttle(feed):
    return next(element for element in feed if isinstance(element, ThrottleValve))


def sample_times(targets, saved_times):
    """Simulation time at which each profile point is judged: the last saved
    time inside its hold (before the next step starts), so start-up and step
    transients are past. None for a point with no saved time in its hold (the
    last point, which shutdown starts from at once, or holds shorter than the
    save interval)."""
    profile_start, profile_end, _ = burn_timeline(targets)
    times = []
    for k, target in enumerate(targets):
        hold_start = profile_start + target["time"]
        hold_end = (profile_start + targets[k + 1]["time"] - step_duration(target, targets[k + 1])
                    if k + 1 < len(targets) else hold_start)
        inside = saved_times[(saved_times > hold_start + 1e-9) & (saved_times <= hold_end + 1e-9)]
        times.append(float(inside[-1]) if inside.size else None)
    return times


def cda_table(filename, targets):
    """Rows (dicts) of the per-profile-point table."""
    time, tracks = LoadTracks(filename)
    at = lambda name, t: float(np.interp(t, time, _track(tracks, name)))
    otv, ftv = _throttle(D.LOX_FEED), _throttle(D.FUEL_FEED)
    rows = []
    for target, t in zip(targets, sample_times(targets, time)):
        if t is None:
            rows.append({"Profile time [s]": round(target["time"], 3), "Sim time [s]": "not held",
                         "Target thrust [N]": round(target["thrust"], 1)})
            continue
        thrust = at("Thrust [N]", t)
        otv_cda, ftv_cda = at(f"{otv.name} CdA [m2]", t), at(f"{ftv.name} CdA [m2]", t)
        rows.append({
            "Profile time [s]": round(target["time"], 3),
            "Sim time [s]": round(t, 3),
            "Target thrust [N]": round(target["thrust"], 1),
            "Thrust [N]": round(thrust, 1),
            "Thrust error [%]": round(100 * (thrust - target["thrust"]) / target["thrust"], 2) if target["thrust"] else 0.0,
            "Pc [psia]": round(at("Chamber Pressure [psia]", t), 1),
            "MR": round(at("Mixture Ratio", t), 3),
            "Ox flow [kg/s]": round(at("Ox Mass Flow [kg/s]", t), 4),
            "Fuel flow [kg/s]": round(at("Fuel Mass Flow [kg/s]", t), 4),
            f"{otv.name} CdA [mm2]": round(otv_cda * 1e6, 2),
            f"{otv.name} % of max": round(100 * otv_cda / otv.max_cda, 1),
            f"{ftv.name} CdA [mm2]": round(ftv_cda * 1e6, 2),
            f"{ftv.name} % of max": round(100 * ftv_cda / ftv.max_cda, 1),
            "LOX injector dP [psi]": round(at("LOX Injector dP [psi]", t), 1),
            "Fuel injector dP [psi]": round(at("Fuel Injector dP [psi]", t), 1),
            "LOX tank [psia]": round(at("LOX Tank Pressure [psia]", t), 1),
            "Fuel tank [psia]": round(at("Fuel Tank Pressure [psia]", t), 1),
            "COPV [psia]": round(at("COPV Pressure [psia]", t), 1),
        })
    return rows


def run_summary(filename, targets, table=None):
    """{metric: value} for one run (also what the Monte Carlo collects)."""
    time, tracks = LoadTracks(filename)
    tr = lambda name: _track(tracks, name)
    table = table if table is not None else cda_table(filename, targets)
    profile_start, profile_end, shutdown_end = burn_timeline(targets)
    in_profile = (time >= profile_start) & (time <= profile_end)

    thrust = tr("Thrust [N]")
    on = thrust >= D.REPORT_THRUST_ON_FRACTION * max(t["thrust"] for t in targets)
    burn_time = float(time[on][-1] - time[on][0]) if on.any() else 0.0
    otv, ftv = _throttle(D.LOX_FEED), _throttle(D.FUEL_FEED)
    pc = tr("Chamber Pressure [psia]")
    judged = [row for row in table if "Thrust [N]" in row]
    errors = np.array([row["Thrust error [%]"] for row in judged])
    relief_names = [b.name for b in D.BRANCHES if "Relief" in b.name]

    summary = {
        "Burn time [s]": round(burn_time, 2),
        "Total impulse [N s]": round(float(np.trapezoid(thrust, time)), 0),
        "LOX residual [kg]": round(float(tr("LOX Liquid Mass [kg]")[-1]), 3),
        "LOX residual [%]": round(100 * float(tr("LOX Liquid Mass [kg]")[-1]) / D.LOX_TANK.liquid_mass, 2),
        "Fuel residual [kg]": round(float(tr("Fuel Liquid Mass [kg]")[-1]), 3),
        "Fuel residual [%]": round(100 * float(tr("Fuel Liquid Mass [kg]")[-1]) / D.FUEL_TANK.liquid_mass, 2),
        "LOX used [kg]": round(D.LOX_TANK.liquid_mass - float(tr("LOX Liquid Mass [kg]")[-1]), 3),
        "Fuel used [kg]": round(D.FUEL_TANK.liquid_mass - float(tr("Fuel Liquid Mass [kg]")[-1]), 3),
        "COPV pressure at shutdown [psia]": round(float(np.interp(shutdown_end, time, tr("COPV Pressure [psia]"))), 1),
        "COPV end pressure [psia]": round(float(tr("COPV Pressure [psia]")[-1]), 1),
        "COPV gas used [kg]": round(float(tr("COPV Gas Mass [kg]")[0] - tr("COPV Gas Mass [kg]")[-1]), 3),
        "COPV gas left [kg]": round(float(tr("COPV Gas Mass [kg]")[-1]), 3),
        "LOX tank min in profile [psia]": round(float(tr("LOX Tank Pressure [psia]")[in_profile].min()), 1),
        "Fuel tank min in profile [psia]": round(float(tr("Fuel Tank Pressure [psia]")[in_profile].min()), 1),
        "LOX tank max [psia]": round(float(tr("LOX Tank Pressure [psia]").max()), 1),
        "Fuel tank max [psia]": round(float(tr("Fuel Tank Pressure [psia]").max()), 1),
        "Min LOX injector stiffness [% Pc]": round(float(np.min(
            tr("LOX Injector dP [psi]")[in_profile] / pc[in_profile]) * 100), 1),
        "Min fuel injector stiffness [% Pc]": round(float(np.min(
            tr("Fuel Injector dP [psi]")[in_profile] / pc[in_profile]) * 100), 1),
        "Max |thrust error| at profile points [%]": round(float(np.max(np.abs(errors))), 2),
        "RMS thrust error at profile points [%]": round(float(np.sqrt(np.mean(errors**2))), 2),
        "Profile points judged": f"{len(judged)} of {len(table)}",
        f"Max {otv.name} % of max CdA": max(row[f"{otv.name} % of max"] for row in judged),
        f"Max {ftv.name} % of max CdA": max(row[f"{ftv.name} % of max"] for row in judged),
        "Relief mass vented [kg]": round(sum(float(np.trapezoid(tr(f"{name} Mass Flow [kg/s]"), time))
                                             for name in relief_names), 4),
    }
    reachable = (summary[f"Max {otv.name} % of max CdA"] < 99.9 and summary[f"Max {ftv.name} % of max CdA"] < 99.9)
    summary["All thrust points reachable"] = "yes" if reachable else "NO (a throttle valve is fully open)"
    summary["Thrust within tolerance"] = ("yes" if summary["Max |thrust error| at profile points [%]"]
                                          <= 100 * D.REPORT_THRUST_TOLERANCE else "NO")
    return summary


def write_report(filename, targets, directory=D.PROP_RESULTS_DIRECTORY, vents=True):
    """Write cda_table.csv, run_summary.csv (and vent_times.csv); returns the summary."""
    os.makedirs(directory, exist_ok=True)
    table = cda_table(filename, targets)
    with open(os.path.join(directory, "cda_table.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(max(table, key=len)), restval="")
        writer.writeheader()
        writer.writerows(table)
    summary = run_summary(filename, targets, table)
    if vents:
        from isolated_sims.Hopper_Vent import vent_times
        for row in vent_times(directory)[0]:
            summary[f"{row['Vent']} time, {row['Case']} tank [s]"] = row["Vent time [s]"]
        from isolated_sims.Hopper_Drain import drain_times
        for row in drain_times(directory)[0]:
            summary[f"LOX drain time at {row['Ullage [psig]']:g} psig [s]"] = row["Drain time [s]"]
        from isolated_sims.Hopper_Tank_Pressure_Sweep import tank_pressure_sweep
        summary.update(tank_pressure_sweep(targets, directory))
        from isolated_sims.Hopper_Purge import purge_check
        summary.update(purge_check(filename, directory))
        from isolated_sims.Hopper_Purge_Sizing import purge_sizing
        summary.update(purge_sizing(filename, directory))
        from isolated_sims.Hopper_Regulator_Sizing import regulator_sizing
        summary.update(regulator_sizing(filename, directory))
    with open(os.path.join(directory, "run_summary.csv"), "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Metric", "Value"])
        writer.writerows(summary.items())
    return summary


if __name__ == "__main__":
    from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
    results = sys.argv[1] if len(sys.argv) > 1 else D.PROP_FILENAME
    summary = write_report(results, HopperThrustTargets(D.PROP_THRUST_CURVE))
    width = max(map(len, summary))
    for key, value in summary.items():
        print(f"{key:<{width}}  {value}")

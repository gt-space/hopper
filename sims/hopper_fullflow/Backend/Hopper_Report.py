"""Single-run report: the numbers asked for after every Hopper run.

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
Vent times come from Hopper_Vent.py (vent_times.csv), LOX drain times from
Hopper_Drain.py (lox_drain_times.csv).

    python Backend/Hopper_Report.py [results.h5]
"""

import csv
import os
import sys

import numpy as np
from scipy.constants import psi

from pathlib import Path

if not __package__:   # run as a script (python Backend/Hopper_Report.py): make HOPPER_DATA and the Backend
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # package importable

import HOPPER_DATA as D
from Backend.Hopper_Elements import ThrottleValve
from Backend.Hopper_Export import LoadTracks, _track
from Backend.Hopper_Prop_System import burn_timeline


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
        hold_end = profile_start + targets[k + 1]["time"] - D.THROTTLE_STEP_TIME if k + 1 < len(targets) else hold_start
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
    profile_start, profile_end, _ = burn_timeline(targets)
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

"""Run metadata and reports: residuals, demonstrated thrust, operating limits."""

import csv
import json
from pathlib import Path

import h5py
import numpy as np

from Hopper_Export import LoadTracks, _track
from Hopper_Inputs import LBF


def write_run_metadata(filename, targets, integration_dt):
    import HOPPER_DATA as data
    profile_start = data.MAIN_VALVE_OPEN_TIME + data.MAIN_VALVE_RAMP_TIME
    metadata = {
        "schema": 1, "profile": data.PROP_THRUST_CURVE_NAME,
        "profile_start": profile_start, "profile_end": profile_start + targets[-1]["time"],
        "integration_dt": integration_dt,
        "hold_time": data.MAX_THRUST_HOLD_TIME,
        "tracking_tolerance": data.MAX_THRUST_TRACKING_TOLERANCE,
        "map_pc_max_pa": data.CHAMBER_MAP_PRESSURE_UPPER_LIMIT,
        "map_mr_min": data.CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT,
        "map_mr_max": data.CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT,
        "copvs": [{"name": b.name, "minimum_temperature": b.min_temperature} for b in data.COPVS],
        "tanks": {"LOX": {"initial_mass": data.LOX_TANK.liquid_mass},
                  "Fuel": {"initial_mass": data.FUEL_TANK.liquid_mass}},
        "throttles": {label: [{"name": e.name, "max_cda": e.max_cda}
                                for e in feed if hasattr(e, "max_cda")]
                      for label, feed in [("LOX", data.LOX_FEED), ("Fuel", data.FUEL_FEED)]},
        "paths": {"ox": [e.name for e in data.LOX_FEED], "fuel": [e.name for e in data.FUEL_FEED],
                  "lox_pressurant": [e.name for e in data.LOX_PRESSURANT_LEG],
                  "fuel_pressurant": [e.name for e in data.FUEL_PRESSURANT_LEG]},
    }
    with h5py.File(filename, "a") as h5:
        h5.attrs["hopper_metadata"] = json.dumps(metadata)


def read_metadata(filename):
    with h5py.File(filename, "r") as h5:
        raw = h5.attrs.get("hopper_metadata")
    return json.loads(raw) if raw else None


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

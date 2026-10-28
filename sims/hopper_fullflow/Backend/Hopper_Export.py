"""CSV exports for the Hopper FullFlow model (plots are in Hopper_Plots.py).

FullFlow writes the transient to HDF5. Each Network.track(...) in
Hopper_Prop_System.py becomes a dataset under
"<network name>/transient/runs/base/tracks", named after the track with spaces
and slashes replaced by underscores. The time vector is in ".../runs/base/time".
The reading helpers follow Vespula_Export.py.
"""

import csv
import re

import h5py
import numpy as np
from scipy.constants import psi


def _safe_group_name(name):
    """Track name -> HDF5 dataset name (same rule FullFlow uses)."""
    name = str(name).strip()
    name = re.sub(r"[\\/\x00]", "_", name)
    name = re.sub(r"\s+", "_", name)
    return name


def _find_network(h5):
    """Name of the first network group in the file that has a transient run."""
    for name, group in h5.items():
        if isinstance(group, h5py.Group) and "transient" in group:
            return name

    raise ValueError("No FullFlow transient network found in HDF5 file.")


def LoadTracks(filename):
    """Return (time, {dataset name: array}) from a FullFlow transient HDF5 file."""
    with h5py.File(filename, "r") as h5:
        run = h5[f"{_find_network(h5)}/transient/runs/base"]
        time = np.asarray(run["time"], dtype=float)
        tracks = {key: np.asarray(run["tracks"][key], dtype=float) for key in run["tracks"]}
    return time, tracks


def _track(tracks, name):
    """Look up a track by the name it was given in Network.track()."""
    key = _safe_group_name(name)
    if key not in tracks:
        raise KeyError(f"Track '{name}' (dataset '{key}') was not found.")
    return tracks[key]


def ExportEngineTargets(targets, filename):
    """Write the targeter output in vehicle_sim.py's engine_targets.csv format."""
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Time [s]", "Target Thrust [N]", "Target Pc [psi]", "Mdot Ox [kg/s]",
                         "Mdot Fu [kg/s]", "Mdot Total [kg/s]"])
        for target in targets:
            total = target["ox_mass_flow"] + target["fuel_mass_flow"]
            writer.writerow([
                target["time"],
                target["thrust"],
                round(target["chamber_pressure"] / psi, 2),
                round(target["ox_mass_flow"], 4),
                round(target["fuel_mass_flow"], 4),
                round(total, 4),
            ])


def ExportActualPerformance(filename, csv_filename):
    """Write the transient engine results in vehicle_sim.py's
    actual_performance.csv format (plus a target-thrust column)."""
    time, tracks = LoadTracks(filename)

    # (CSV header, track name, decimals)
    columns = [
        ("Actual Pc [psi]", "Chamber Pressure [psia]", 2),
        ("Actual Thrust [N]", "Thrust [N]", 2),
        ("Target Thrust [N]", "Target Thrust [N]", 2),
        ("Actual MR", "Mixture Ratio", 4),
        ("Mdot Ox [kg/s]", "Ox Mass Flow [kg/s]", 4),
        ("Mdot Fu [kg/s]", "Fuel Mass Flow [kg/s]", 4),
        ("Total Mdot [kg/s]", "Total Mass Flow [kg/s]", 4),
        ("Actual C* [m/s]", "C* [m/s]", 2),
        ("Actual Isp [s]", "Isp [s]", 2),
    ]
    data = [_track(tracks, track) for _, track, _ in columns]

    # C* and Isp are flow ratios, meaningless before ignition and after tail-off.
    # As in vehicle_sim.Engine, report zero unless Pc > 1.5 x ambient.
    pc = _track(tracks, "Chamber Pressure [psia]")
    lit = pc > 1.5 * pc[0]
    for i, (header, _, _) in enumerate(columns):
        if header in ("Actual C* [m/s]", "Actual Isp [s]"):
            data[i] = np.where(lit, data[i], 0.0)

    with open(csv_filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Time [s]"] + [header for header, _, _ in columns])
        for i, t in enumerate(time):
            writer.writerow([round(t, 3)] + [round(values[i], digits)
                                            for values, (_, _, digits) in zip(data, columns)])

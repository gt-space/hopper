"""Readers for input files referenced from HOPPER_DATA.py.

load_thrust_curve(path)
    A thrust curve exported as CSV (e.g. from the trajectory sim). The file
    needs a header row; the time column is the first whose header contains
    "time" (or "t"), the thrust column the first containing "thrust" (or
    "force"/"f"). Units are read from the header: time in s (or ms), thrust in
    N (or kN, lbf). Other columns are ignored. Rows that are not numbers are
    skipped. Time is shifted so the curve starts at 0.

    Trajectory exports are often dense (every 0.01 s or finer). Every curve
    point becomes a breakpoint the transient solver must land on, so a curve
    with points closer than min_spacing is resampled onto a grid of that
    spacing (linear interpolation); run such curves with
    PROP_PROFILE_INTERPOLATION = "linear".
"""

import csv
from pathlib import Path

import numpy as np
from scipy.constants import kilo, milli, pound_force

PROJECT = Path(__file__).resolve().parents[1]   # hopper_fullflow; curve paths are relative to it


def _column(headers, keys):
    for i, header in enumerate(headers):
        name = header.strip().lower()
        if any(name == key or name.startswith(key + " ") or name.startswith(key + "[") or key in name.split()
               for key in keys):
            return i
    for i, header in enumerate(headers):
        if any(key in header.strip().lower() for key in keys if len(key) > 1):
            return i
    raise ValueError(f"No column for {keys} in header {headers}")


def load_thrust_curve(path, min_spacing=0.1):
    """{profile time [s]: thrust [N]} from a CSV file (path relative to the hopper_fullflow folder)."""
    path = Path(path)
    if not path.is_absolute():
        path = PROJECT / path
    with open(path, newline="", encoding="utf-8-sig") as file:
        rows = list(csv.reader(file))
    headers = rows[0]
    t_col = _column(headers, ["time", "t"])
    f_col = _column(headers, ["thrust", "force", "f"])
    t_header, f_header = headers[t_col].lower(), headers[f_col].lower()
    t_scale = milli if "ms" in t_header else 1.0
    f_scale = pound_force if "lbf" in f_header else kilo if "kn" in f_header else 1.0

    times, thrusts = [], []
    for row in rows[1:]:
        try:
            times.append(float(row[t_col]) * t_scale)
            thrusts.append(float(row[f_col]) * f_scale)
        except (ValueError, IndexError):
            continue
    if len(times) < 2:
        raise ValueError(f"{path}: fewer than two numeric rows")
    order = np.argsort(times)
    times = np.asarray(times)[order] - min(times)
    thrusts = np.asarray(thrusts)[order]

    if np.min(np.diff(times)) < min_spacing - 1e-9:
        grid = np.arange(0.0, times[-1] + 1e-9, min_spacing)
        if grid[-1] < times[-1] - 1e-9:
            grid = np.append(grid, times[-1])
        thrusts = np.interp(grid, times, thrusts)
        times = grid
    return {round(float(t), 6): float(f) for t, f in zip(times, thrusts)}

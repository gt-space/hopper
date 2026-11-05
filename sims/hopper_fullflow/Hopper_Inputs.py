"""Validated numerical inputs shared by configuration, physics, and tests."""

import csv
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
LBF = 4.4482216152605


def ramp_profile(start_thrust, end_thrust, duration, points=101):
    """Linearly sampled thrust demand, N and s; both ramp directions are allowed."""
    if not all(math.isfinite(float(v)) for v in (start_thrust, end_thrust, duration)):
        raise ValueError("Ramp inputs must be finite.")
    if min(start_thrust, end_thrust) <= 0 or duration <= 0 or int(points) != points or points < 2:
        raise ValueError("Ramp requires positive thrust/duration and at least two points.")
    return dict(zip(np.linspace(0, duration, int(points)).tolist(),
                    np.linspace(start_thrust, end_thrust, int(points)).tolist()))


def validate_profile(profile):
    points = sorted((float(t), float(f)) for t, f in profile.items())
    if not points or points[0][0] != 0 or any(not math.isfinite(t + f) or f <= 0 for t, f in points):
        raise ValueError("Thrust profile must start at t=0 with finite, positive thrusts.")
    if any(b[0] <= a[0] for a, b in zip(points, points[1:])):
        raise ValueError("Thrust-profile times must be distinct and increasing.")
    return dict(points)


@lru_cache(maxsize=32)
def _read_table(filename, modified_ns, value_column):
    with open(filename, newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        headers = reader.fieldnames or []
        if "pc_psia" not in headers or value_column not in headers:
            raise ValueError(f"{filename}: requires pc_psia and {value_column} columns.")
        has_mr = "mr" in headers
        rows = [(float(row["pc_psia"]), float(row["mr"]) if has_mr else 0.0,
                 float(row[value_column])) for row in reader]
    if not rows or not all(math.isfinite(x) for row in rows for x in row):
        raise ValueError(f"{filename}: table must contain finite numeric data.")
    if any(p < 0 or (has_mr and m <= 0) for p, m, _ in rows):
        raise ValueError(f"{filename}: Pc must be nonnegative and MR positive.")
    if value_column == "efficiency" and any(not 0 < value <= 1 for _, _, value in rows):
        raise ValueError(f"{filename}: every efficiency must be in (0, 1].")
    if len(set((p, m) for p, m, _ in rows)) != len(rows):
        raise ValueError(f"{filename}: duplicate Pc/MR points.")
    pressures = sorted(set(p for p, _, _ in rows))
    ratios = sorted(set(m for _, m, _ in rows))
    lookup = {(p, m): value for p, m, value in rows}
    if len(rows) != len(pressures) * len(ratios):
        raise ValueError(f"{filename}: Pc/MR data must form a complete rectangular grid.")
    values = np.array([[lookup[p, m] for p in pressures] for m in ratios])
    return pressures, ratios, values


def table_value(filename, pc_pa, mr, value_column="efficiency"):
    """Linear/bilinear interpolation; clamp inputs at table boundaries."""
    path = Path(filename)
    if not path.is_absolute():
        path = ROOT / path
    pressures, ratios, values = _read_table(str(path), path.stat().st_mtime_ns, value_column)
    along_pc = [np.interp(pc_pa / 6894.757293168, pressures, row) for row in values]
    return float(np.interp(mr, ratios, along_pc))


def efficiency_value(function, filename, pc_pa, mr):
    value = table_value(filename, pc_pa, mr) if filename else float(function(pc_pa, mr))
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError(f"Engine efficiency must be in (0, 1], got {value} at Pc={pc_pa}, MR={mr}.")
    return value


def delta_t_value(points, pc_pa):
    points = sorted((float(p), float(dt)) for p, dt in points)
    if not points or any(not math.isfinite(p + dt) or p < 0 or dt < 0 for p, dt in points):
        raise ValueError("Regen curve requires finite, nonnegative (Pc [Pa], dT [K]) points.")
    if any(b[0] <= a[0] for a, b in zip(points, points[1:])):
        raise ValueError("Regen Pc points must be distinct.")
    return float(np.interp(pc_pa, [p for p, _ in points], [dt for _, dt in points]))


def slew_position(previous, command, elapsed, open_time, close_time):
    """Full-stroke linear travel, with different opening and closing rates."""
    previous = min(1.0, max(0.0, float(previous)))
    command = min(1.0, max(0.0, float(command)))
    travel_time = open_time if command >= previous else close_time
    if travel_time == 0:
        return command
    change = max(0.0, elapsed) / travel_time
    return previous + min(change, max(-change, command - previous))

"""Compare the old general_fluid_network vehicle sim with the FullFlow model.

Inputs
  --old-history  vehicle_sim_history.csv written by run_vehicle_sim.py
  --new-results  Hopper_Prop_System.h5 written by Hopper.py
  --output-dir   where the comparison plots and summary go

Time alignment
  Both models fly the same thrust profile, but on different clocks:
    * vehicle_sim.py starts the profile at t = 0 with the engine already at
      full flow (no start-up), and runs to t = 29.6 s.
    * The FullFlow model ramps the main valves open first; its
      profile time 0 is at MAIN_VALVE_OPEN_TIME + MAIN_VALVE_RAMP_TIME (0.5 s).
  Everything here is plotted against profile time: the new model's time minus
  0.5 s. The difference table uses profile time 1.0-29.5 s, which skips
  vehicle_sim.py's first second (its chamber pressure relaxes up from 67 psia)
  and the FullFlow start-up.

Outputs
  comparison_*.png         old vs new overlays
  comparison_summary.csv   mean / max differences for each quantity
  comparison_summary.md    the same table, plus end-of-burn values
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# The FullFlow model's modules live one folder up.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from HOPPER_DATA import MAIN_VALVE_OPEN_TIME, MAIN_VALVE_RAMP_TIME, THRUST_CURVE_1DOF
from Hopper_Export import LoadTracks, _track


PROFILE_START = MAIN_VALVE_OPEN_TIME + MAIN_VALVE_RAMP_TIME   # new-model time of profile t = 0
WINDOW = (1.0, 29.5)                                          # profile-time window for the statistics
OLD_COLOR, NEW_COLOR, TARGET_COLOR = "tab:blue", "tab:orange", "black"

# Quantities compared: (label, column/track name shared by both files, units).
# run_vehicle_sim.py names its CSV columns after the FullFlow tracks, so one name
# addresses both.
QUANTITIES = [
    ("Thrust", "Thrust [N]", "N"),
    ("Chamber pressure", "Chamber Pressure [psia]", "psia"),
    ("Mixture ratio", "Mixture Ratio", "-"),
    ("Ox mass flow", "Ox Mass Flow [kg/s]", "kg/s"),
    ("Fuel mass flow", "Fuel Mass Flow [kg/s]", "kg/s"),
    ("Isp", "Isp [s]", "s"),
    ("C*", "C* [m/s]", "m/s"),
    ("LOX tank pressure", "LOX Tank Pressure [psia]", "psia"),
    ("Fuel tank pressure", "Fuel Tank Pressure [psia]", "psia"),
    ("COPV pressure", "COPV Pressure [psia]", "psia"),
    ("COPV temperature", "COPV Temperature [K]", "K"),
    ("LOX ullage temperature", "LOX Ullage Temperature [K]", "K"),
    ("Fuel ullage temperature", "Fuel Ullage Temperature [K]", "K"),
    ("LOX injector dP", "LOX Injector dP [psi]", "psi"),
    ("Fuel injector dP", "Fuel Injector dP [psi]", "psi"),
    ("Regen dP", "Regen Jacket dP [psi]", "psi"),
    ("OTV CdA", "OTV CdA [mm2]", "mm^2"),
    ("FTV CdA", "FTV CdA [mm2]", "mm^2"),
    ("LOX remaining", "LOX Liquid Mass [kg]", "kg"),
    ("Fuel remaining", "Fuel Liquid Mass [kg]", "kg"),
    ("COPV gas remaining", "COPV Gas Mass [kg]", "kg"),
]


def load_old(filename):
    """vehicle_sim_history.csv -> (time, {column: array}). Time is already profile time."""
    with open(filename, newline="") as file:
        rows = list(csv.reader(file))
    header, data = rows[0], np.array(rows[1:], dtype=float)
    columns = {name: data[:, i] for i, name in enumerate(header)}
    return columns["Time [s]"], columns


def load_new(filename):
    """FullFlow HDF5 -> (profile time, {track name: array})."""
    time, tracks = LoadTracks(filename)
    columns = {name: _track(tracks, name) for _, name, _ in QUANTITIES}
    columns["Target Thrust [N]"] = _track(tracks, "Target Thrust [N]")
    return time - PROFILE_START, columns


def staircase_target(time):
    """Thrust commanded at each profile time: the most recent profile point.

    Both models command their valves this way (vehicle_sim.py's action dict,
    and the FullFlow model with PROP_PROFILE_INTERPOLATION = "step").
    """
    points = np.array(list(THRUST_CURVE_1DOF.keys()))
    values = np.array(list(THRUST_CURVE_1DOF.values()))
    # Round first: a sample saved at a step time can read 5.4999999 s, which
    # would otherwise pick the previous point.
    index = np.clip(np.searchsorted(points, np.round(time, 6), side="right") - 1, 0, len(points) - 1)
    return values[index]


def plateau_times():
    """The last 0.1 s sample before each profile point, inside WINDOW.

    Each profile point starts a hold. Right after a step the two models differ
    mostly because vehicle_sim.py relaxes Pc 30 % per 0.1 s step toward its new
    value while the FullFlow chamber responds at once. The end of each hold is
    where both have settled the most, so it isolates the steady model physics.
    """
    points = list(THRUST_CURVE_1DOF.keys())
    times = [round(next_point - 0.1, 1) for point, next_point in zip(points[:-1], points[1:])
             if next_point - 0.1 > point + 1e-9]
    return np.array([t for t in times if WINDOW[0] <= t <= WINDOW[1]])


def statistics(old_time, old, new_time, new):
    """Per-quantity differences on the old model's 0.1 s grid inside WINDOW,
    over every sample and over the hold plateaus only."""
    mask = (old_time >= WINDOW[0]) & (old_time <= WINDOW[1])
    grid = old_time[mask]
    plateaus = plateau_times()
    rows = []
    for label, name, units in QUANTITIES:
        old_values = old[name][mask]
        new_values = np.interp(grid, new_time, new[name])
        difference = new_values - old_values
        scale = np.maximum(np.abs(old_values), 1e-12)

        old_plateau = np.interp(plateaus, old_time, old[name])
        new_plateau = np.interp(plateaus, new_time, new[name])
        plateau_error = np.abs(new_plateau - old_plateau) / np.maximum(np.abs(old_plateau), 1e-12)

        rows.append({
            "Quantity": label,
            "Units": units,
            "Old mean": np.mean(old_values),
            "New mean": np.mean(new_values),
            "Mean difference (new - old)": np.mean(difference),
            "Mean |difference| [%]": 100 * np.mean(np.abs(difference) / scale),
            "Max |difference| [%]": 100 * np.max(np.abs(difference) / scale),
            "Plateau mean |difference| [%]": 100 * np.mean(plateau_error),
            "Plateau max |difference| [%]": 100 * np.max(plateau_error),
        })
    return rows


def thrust_tracking(time, thrust):
    """RMS and max |thrust - commanded thrust| inside WINDOW [N]."""
    mask = (time >= WINDOW[0]) & (time <= WINDOW[1])
    error = thrust[mask] - staircase_target(time[mask])
    return float(np.sqrt(np.mean(error**2))), float(np.max(np.abs(error)))


def end_values(time, data, at=WINDOW[1]):
    """Value of each quantity at profile time `at`."""
    return {name: float(np.interp(at, time, data[name])) for _, name, _ in QUANTITIES}


def overlay(ax, old_time, old_values, new_time, new_values, ylabel, target=None):
    ax.plot(old_time, old_values, color=OLD_COLOR, label="vehicle_sim.py (general_fluid_network)")
    ax.plot(new_time, new_values, color=NEW_COLOR, label="FullFlow model")
    if target is not None:
        ax.plot(*target, color=TARGET_COLOR, linestyle="--", linewidth=1, label="Profile target")
    ax.set_ylabel(ylabel)
    ax.grid(True)


def plot_comparison(old_time, old, new_time, new, output_dir):
    """One PNG per group of related quantities."""
    groups = [
        ("comparison_engine.png", "Engine: vehicle_sim.py vs FullFlow",
         ["Thrust [N]", "Chamber Pressure [psia]", "Mixture Ratio"]),
        ("comparison_flows.png", "Propellant flows and throttle-valve CdA",
         ["Ox Mass Flow [kg/s]", "Fuel Mass Flow [kg/s]", "OTV CdA [mm2]", "FTV CdA [mm2]"]),
        ("comparison_pressurization.png", "Tanks and COPV",
         ["LOX Tank Pressure [psia]", "Fuel Tank Pressure [psia]", "COPV Pressure [psia]"]),
        ("comparison_temperatures.png", "Gas temperatures",
         ["COPV Temperature [K]", "LOX Ullage Temperature [K]", "Fuel Ullage Temperature [K]"]),
        ("comparison_injectors.png", "Injector and regen pressure drops",
         ["LOX Injector dP [psi]", "Fuel Injector dP [psi]", "Regen Jacket dP [psi]"]),
        ("comparison_propellant.png", "Remaining propellant and pressurant",
         ["LOX Liquid Mass [kg]", "Fuel Liquid Mass [kg]", "COPV Gas Mass [kg]"]),
    ]
    profile_time = np.array(list(THRUST_CURVE_1DOF.keys()))
    profile_thrust = np.array(list(THRUST_CURVE_1DOF.values()))

    for filename, title, names in groups:
        fig, axs = plt.subplots(len(names), 1, figsize=(11, 2.6 * len(names) + 1), sharex=True)
        for ax, name in zip(np.atleast_1d(axs), names):
            target = (profile_time, profile_thrust) if name == "Thrust [N]" else None
            overlay(ax, old_time, old[name], new_time, new[name], name, target)
        np.atleast_1d(axs)[0].legend(fontsize=8, loc="best")
        np.atleast_1d(axs)[-1].set_xlabel("Profile time [s]  (FullFlow time - 0.5 s)")
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, filename), dpi=150)
        plt.close(fig)

    # Injector stiffness (dP / Pc), the figure vehicle_sim.py makes, for both models.
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for time, data, color, label in [(old_time, old, OLD_COLOR, "vehicle_sim.py"),
                                     (new_time, new, NEW_COLOR, "FullFlow")]:
        lit = data["Chamber Pressure [psia]"] > 1.5 * 14.8
        pc = data["Chamber Pressure [psia]"][lit]
        ax.plot(time[lit], 100 * data["LOX Injector dP [psi]"][lit] / pc, color=color, label=f"OINJ, {label}")
        ax.plot(time[lit], 100 * data["Fuel Injector dP [psi]"][lit] / pc, color=color, linestyle=":",
                label=f"FINJ, {label}")
    ax.set_xlabel("Profile time [s]")
    ax.set_ylabel("Stiffness [%]")
    ax.set_title("Injector stiffness")
    ax.grid(True)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "comparison_stiffness.png"), dpi=150)
    plt.close(fig)


def write_summary(rows, tracking, ends, output_dir):
    """comparison_summary.csv and comparison_summary.md."""
    keys = list(rows[0].keys())
    with open(os.path.join(output_dir, "comparison_summary.csv"), "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

    def fmt(value):
        return f"{value:.4g}" if isinstance(value, float) else str(value)

    lines = [
        "# vehicle_sim.py vs FullFlow model (1dof thrust profile)",
        "",
        f"Statistics over profile time {WINDOW[0]}-{WINDOW[1]} s, on vehicle_sim.py's 0.1 s grid.",
        "\"Plateau\" columns use only the last sample of each profile hold "
        f"({len(plateau_times())} points), where both models have settled after a step.",
        "",
        "| " + " | ".join(keys) + " |",
        "|" + "---|" * len(keys),
    ]
    lines += ["| " + " | ".join(fmt(row[key]) for key in keys) + " |" for row in rows]
    lines += [
        "",
        "## Thrust vs commanded thrust (most recent profile point)",
        "",
        "| Model | RMS error [N] | Max error [N] |",
        "|---|---|---|",
    ]
    lines += [f"| {model} | {rms:.1f} | {peak:.1f} |" for model, (rms, peak) in tracking.items()]
    lines += [
        "",
        f"## Values at profile time {WINDOW[1]} s",
        "",
        "| Quantity | vehicle_sim.py | FullFlow |",
        "|---|---|---|",
    ]
    lines += [f"| {label} [{units}] | {ends['old'][name]:.4g} | {ends['new'][name]:.4g} |"
              for label, name, units in QUANTITIES]
    with open(os.path.join(output_dir, "comparison_summary.md"), "w") as file:
        file.write("\n".join(lines) + "\n")
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--old-history", required=True)
    parser.add_argument("--new-results", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    old_time, old = load_old(args.old_history)
    new_time, new = load_new(args.new_results)

    plot_comparison(old_time, old, new_time, new, args.output_dir)
    rows = statistics(old_time, old, new_time, new)
    tracking = {"vehicle_sim.py": thrust_tracking(old_time, old["Thrust [N]"]),
                "FullFlow": thrust_tracking(new_time, new["Thrust [N]"])}
    ends = {"old": end_values(old_time, old), "new": end_values(new_time, new)}
    print(write_summary(rows, tracking, ends, args.output_dir))

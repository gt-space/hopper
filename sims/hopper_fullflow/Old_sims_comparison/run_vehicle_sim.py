"""Run the old general_fluid_network vehicle sim and save its outputs.

This runs hopper/sims/general_fluid_network/vehicle_sim.py exactly as it is in
the repository, with two in-memory changes (the file on disk is never edited):

  1. Thrust profile. The script selects
         thrust_curve = thrust_curve_throttle_test
     and that line is swapped for
         thrust_curve = thrust_curve_1dof
     so both models fly the full 29.7 s 1dof profile.

  2. Plotting. matplotlib is switched to the non-interactive "Agg" backend
     before the script runs, so its plt.show() calls return immediately
     instead of opening windows. Every figure the script makes is then
     saved as a PNG.

Everything else (hardware values, targeter, efficiency curves, dt = 0.1 s,
explicit integration, CSV exports) is the original code.

Outputs, written to --output-dir:
  engine_targets.csv          the script's own export (unchanged)
  actual_performance.csv      the script's own export (unchanged)
  vehicle_sim_history.csv     every history the script's objects keep
                              (tanks, COPV, valves, injectors, engine), one
                              row per 0.1 s step, for the comparison
  *.png                       the script's own figures

Usage (from any folder, with the general_fluid_network requirements installed):
  python run_vehicle_sim.py --gfn-dir <path to sims/general_fluid_network> --output-dir <folder>
"""

import argparse
import csv
import os
import re
import sys
import time

# Must happen before vehicle_sim.py imports pyplot, otherwise plt.show() blocks.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PSI_TO_PA = 6894.75729   # same conversion vehicle_sim.py uses

# The one line of vehicle_sim.py that is changed, and what it becomes.
PROFILE_LINE = "thrust_curve = thrust_curve_throttle_test"
PROFILE_REPLACEMENT = "thrust_curve = thrust_curve_1dof"


def run_vehicle_sim(gfn_dir, output_dir):
    """Execute vehicle_sim.py with the profile swap and return its globals."""
    script_path = os.path.join(gfn_dir, "vehicle_sim.py")
    with open(script_path, encoding="utf-8") as file:
        source = file.read()

    # Fail loudly if the script has changed and the swap no longer applies,
    # rather than silently running the wrong profile.
    if source.count(PROFILE_LINE) != 1:
        raise RuntimeError(f"Expected exactly one '{PROFILE_LINE}' line in {script_path}.")
    source = source.replace(PROFILE_LINE, PROFILE_REPLACEMENT)

    # vehicle_sim.py imports general_fluid_network from its own folder and writes
    # its CSVs to the current directory, so run it from the output folder with
    # the gfn folder on the import path.
    sys.path.insert(0, os.path.abspath(gfn_dir))
    os.makedirs(output_dir, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(output_dir)
    try:
        namespace = {"__name__": "__main__", "__file__": script_path}
        exec(compile(source, script_path, "exec"), namespace)
    finally:
        os.chdir(cwd)
    return namespace


def save_figures(output_dir):
    """Save every open matplotlib figure, named after its title."""
    for number in plt.get_fignums():
        figure = plt.figure(number)
        title = figure._suptitle.get_text() if figure._suptitle else ""
        if not title and figure.axes:
            title = figure.axes[0].get_title()
        slug = re.sub(r"[^a-z0-9]+", "_", (title or f"figure_{number}").lower()).strip("_")
        figure.savefig(os.path.join(output_dir, f"vehicle_sim_{number:02d}_{slug}.png"), dpi=150)


def export_history(namespace, filename):
    """Write the per-step histories kept by vehicle_sim.py's objects to one CSV.

    general_fluid_network logs every node and connection once per step at the
    same time stamps, so the columns line up row for row.
    """
    engine = namespace["engine"]
    copv = namespace["copv"]
    ox_tank, fu_tank = namespace["ox_tank"], namespace["fu_tank"]
    obb, fbb = namespace["obb"], namespace["fbb"]
    obb_series, fbb_series = namespace["obb_series"], namespace["fbb_series"]
    otv, ftv = namespace["otv"], namespace["ftv"]
    otv_series, ftv_series = namespace["otv_series"], namespace["ftv_series"]
    ox_inj, fu_inj, fu_regen = namespace["ox_inj"], namespace["fu_inj"], namespace["fu_regen"]

    # (column name, history array). Pressures in psia, flows in kg/s, CdA in mm^2.
    columns = [
        ("Time [s]", engine.history["time"]),
        ("Chamber Pressure [psia]", np.array(engine.history["P"]) / PSI_TO_PA),
        ("Thrust [N]", engine.history["thrust"]),
        ("Mixture Ratio", engine.history["MR"]),
        ("Ox Mass Flow [kg/s]", engine.history["mdot_ox"]),
        ("Fuel Mass Flow [kg/s]", engine.history["mdot_fu"]),
        ("C* [m/s]", engine.history["cstar"]),
        ("Isp [s]", engine.history["Isp"]),
        ("Chamber Temperature [K]", engine.history["T"]),
        ("COPV Pressure [psia]", np.array(copv.history["P"]) / PSI_TO_PA),
        ("COPV Temperature [K]", copv.history["T"]),
        ("COPV Gas Mass [kg]", copv.history["m"]),
        ("LOX Tank Pressure [psia]", np.array(ox_tank.history["P"]) / PSI_TO_PA),
        ("LOX Liquid Mass [kg]", ox_tank.history["m"]),
        ("LOX Ullage Temperature [K]", ox_tank.ullage.history["T"]),
        ("LOX Ullage Gas Mass [kg]", ox_tank.ullage.history["m"]),
        ("Fuel Tank Pressure [psia]", np.array(fu_tank.history["P"]) / PSI_TO_PA),
        ("Fuel Liquid Mass [kg]", fu_tank.history["m"]),
        ("Fuel Ullage Temperature [K]", fu_tank.ullage.history["T"]),
        ("Fuel Ullage Gas Mass [kg]", fu_tank.ullage.history["m"]),
        ("OBB State", np.array(obb.history["state"], dtype=float)),
        ("FBB State", np.array(fbb.history["state"], dtype=float)),
        ("OBB Mass Flow [kg/s]", obb_series.history["mdot"]),
        ("FBB Mass Flow [kg/s]", fbb_series.history["mdot"]),
        ("OTV CdA [mm2]", np.array(otv.history["CdA"]) * 1e6),
        ("FTV CdA [mm2]", np.array(ftv.history["CdA"]) * 1e6),
        ("OTV Series Mass Flow [kg/s]", otv_series.history["mdot"]),
        ("FTV Series Mass Flow [kg/s]", ftv_series.history["mdot"]),
        ("LOX Injector dP [psi]", np.array(ox_inj.history["dP"]) / PSI_TO_PA),
        ("Fuel Injector dP [psi]", np.array(fu_inj.history["dP"]) / PSI_TO_PA),
        ("Regen dP [psi]", np.array(fu_regen.history["dP"]) / PSI_TO_PA),
    ]

    rows = min(len(values) for _, values in columns)
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([name for name, _ in columns])
        for i in range(rows):
            writer.writerow([float(values[i]) for _, values in columns])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gfn-dir", required=True, help="path to hopper/sims/general_fluid_network")
    parser.add_argument("--output-dir", required=True, help="folder for the old model's outputs")
    args = parser.parse_args()

    start = time.perf_counter()
    namespace = run_vehicle_sim(args.gfn_dir, args.output_dir)
    save_figures(args.output_dir)
    export_history(namespace, os.path.join(args.output_dir, "vehicle_sim_history.csv"))
    print(f"vehicle_sim.py (1dof profile) finished in {time.perf_counter() - start:.1f} s; "
          f"outputs in {args.output_dir}")

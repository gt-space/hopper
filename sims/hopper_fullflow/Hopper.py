"""Run the Hopper vehicle fluids model.

Every input is in hopper_config.json; HOPPER_DATA.py loads it into Python names.

How it runs (each step can be switched in the JSON):
  1. Chamber map: build Backend/Hopper_Chamber_Map.h5 if asked to or if it is missing.
  2. Thrust targets: Pc and propellant flows for every thrust-profile point.
  3. Prop system: the FullFlow transient of the whole vehicle P&ID.
  4. Results: CSVs (vehicle_sim.py formats) and the plots listed in outputs.plots.
  5. Report (outputs.report): run_summary.csv (burn time, residuals, COPV use,
     margins, thrust checks), cda_table.csv (OTV/FTV CdA per thrust point),
     component_table.csv/.md, schematic.png and Hopper_Parameters.xlsx.

Each solve is archived with its exact inputs in Hopper_Results/<UTC time>/, and the
results are copied to Outputs/ for viewing.

Run from this folder:  python Hopper.py
    python Hopper.py --config path/to/case.json      another case
    python Hopper.py --results-only                  redo the outputs from Outputs/Hopper_Prop_System.h5
Plots can also be made later from the saved results without re-running:
    python Backend/Hopper_Plots.py list                      (everything that can be plotted)
    python Backend/Hopper_Plots.py plot ox_ladder engine     (named plots)
    python Backend/Hopper_Plots.py nodes "LOX Tank" COPV     (any nodes, elements, tracks...)
"""

import argparse
import json
import os
import time

run_started = time.perf_counter()

import matplotlib

parser = argparse.ArgumentParser(description='Run Hopper from a JSON configuration.')
parser.add_argument('--config', help='JSON input file (default: hopper_config.json or HOPPER_CONFIG).')
parser.add_argument('--results-only', action='store_true', help='Redo the outputs from an existing history, without solving.')
args = parser.parse_args()
if args.config:
    os.environ['HOPPER_CONFIG'] = args.config

import HOPPER_DATA as data

if args.results_only:
    data.PROP_SOLVE_PROP_SYSTEM = False
if data.PROP_SOLVE_PROP_SYSTEM:
    for tank in (data.LOX_TANK, data.FUEL_TANK):
        print(f"{tank.name}: initial {tank.pressure / data.psi:g} psia, "
              f"setpoint {tank.set_pressure / data.psi:g} psia "
              f"(-{tank.lower_limit / data.psi:g}/+{tank.upper_limit / data.psi:g} psi).", flush=True)

from Backend.Hopper_Run_Config import finish_run, prepare_run, write_resolved

if data.PROP_SOLVE_PROP_SYSTEM:
    prepare_run(data)

from HOPPER_DATA import *

if not PROP_SOLVE_PROP_SYSTEM and not os.path.exists(PROP_FILENAME):
    raise SystemExit(f"No saved history at {PROP_FILENAME}. Only the latest run keeps one (in Outputs/), unless "
                     "outputs.keep_run_history is true; rerun without --results-only to make it again.")

from Backend.Hopper_Chamber_Map import HopperChamberMap
from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
from Backend.Hopper_Prop_System import HopperPropSystem
from Backend.Hopper_Export import ExportEngineTargets, ExportActualPerformance
from Backend.Hopper_Report import export_summary, read_metadata
from Backend.Hopper_Diagram import hardware_diagram




# Render plots off-screen unless they are to be shown (must precede any pyplot import).
if not PROP_SHOW_PLOTS:
    matplotlib.use("Agg")

from Backend.Hopper_Plots import HopperRun, make_plots




# ---- Chamber Map ---- #
if PROP_GENERATE_CHAMBER_MAP or not os.path.exists(CHAMBER_MAP_FILENAME):

    print("Generating chamber Pc/MR map...")

    HopperChamberMap(filename=CHAMBER_MAP_FILENAME)

if PROP_SOLVE_PROP_SYSTEM:
    write_resolved(data)




# ---- Thrust Targets ---- #
# The flows OTV and FTV will be commanded to, one set per profile point.
print("Solving thrust targets...")

targets = HopperThrustTargets(
    PROP_THRUST_CURVE,
    chamber_map_filename=CHAMBER_MAP_FILENAME,
)

os.makedirs(PROP_RESULTS_DIRECTORY, exist_ok=True)
if PROP_COMPONENT_DIAGRAM:
    hardware_diagram(PROP_RESULTS_DIRECTORY)
ExportEngineTargets(targets, os.path.join(PROP_RESULTS_DIRECTORY, "engine_targets.csv"))




# ---- Prop System ---- #
if PROP_SOLVE_PROP_SYSTEM:

    print("Solving propulsion system...")
    start = time.perf_counter()

    HopperPropSystem(
        targets,
        dt=PROP_DT,
        save_dt=PROP_SAVE_DT,
        filename=PROP_FILENAME,
        chamber_map_filename=CHAMBER_MAP_FILENAME,
        profile_interpolation=PROP_PROFILE_INTERPOLATION,
        rtol=PROP_RTOL,
        jacobian_method=PROP_JACOBIAN_METHOD,
        verbose=PROP_VERBOSE,
        statistics=PROP_STATISTICS,
    )

    print(f"Propulsion system solved in {time.perf_counter() - start:.1f} s.")

if PROP_SOLVE_PROP_SYSTEM or PROP_EXPORT_RESULTS:
    if read_metadata(PROP_FILENAME) is not None:
        export_summary(PROP_FILENAME, PROP_RESULTS_DIRECTORY)
    else:
        print("Saved run predates summary metadata; rerun to obtain maximum-thrust and residual reports.")




# ---- Results ---- #
if PROP_EXPORT_RESULTS:

    print("Exporting results...")

    ExportActualPerformance(PROP_FILENAME, os.path.join(PROP_RESULTS_DIRECTORY, "actual_performance.csv"))
    make_plots(
        HopperRun(PROP_FILENAME),
        PROP_PLOTS,
        PROP_RESULTS_DIRECTORY,
        units=PROP_PLOT_UNITS,
        csv_out=PROP_PLOT_CSV,
        profile_time=PROP_PLOT_PROFILE_TIME,
        show=PROP_SHOW_PLOTS,
    )

    print(f"Results written to {PROP_RESULTS_DIRECTORY}/")




# ---- Report ---- #
if PROP_REPORT:

    print("Writing report...")

    from Backend.Hopper_Report import write_report
    from Backend.Hopper_Hardware import write_table
    from Backend.Hopper_Schematic import draw
    from Backend.Hopper_Parameters import write_workbook

    summary = write_report(PROP_FILENAME, targets, PROP_RESULTS_DIRECTORY)
    write_table(PROP_RESULTS_DIRECTORY)
    draw(os.path.join(PROP_RESULTS_DIRECTORY, "schematic.png"))
    write_workbook(PROP_RESULTS_DIRECTORY)

    width = max(map(len, summary))
    for key, value in summary.items():
        print(f"  {key:<{width}}  {value}")




# ---- Finish ---- #
if PROP_SOLVE_PROP_SYSTEM:
    elapsed = time.perf_counter() - run_started
    with open(os.path.join(PROP_RESULTS_DIRECTORY, 'run_timing.json'), 'w', encoding='utf-8') as stream:
        json.dump({'total_wall_time_s': elapsed}, stream, indent=2)
    finish_run(data)
    print(f"Complete run including setup, solve, and outputs: {elapsed:.1f} s ({elapsed / 60:.2f} min).")

"""Run the Hopper vehicle fluids model.

Every input, including these run settings, is in HOPPER_DATA.py.

Steps, each switchable in HOPPER_DATA.py:
  1. Chamber map: build Hopper_Chamber_Map.h5 if asked to or if it is missing.
  2. Thrust targets: Pc and propellant flows for every thrust-profile point.
  3. Prop system: the FullFlow transient of the whole vehicle P&ID.
  4. Results: CSVs (vehicle_sim.py formats) and the plots listed in PROP_PLOTS
     in PROP_RESULTS_DIRECTORY.
  5. Report (PROP_REPORT): run_summary.csv (burn time, residuals, COPV use,
     margins, thrust checks), cda_table.csv (OTV/FTV CdA per thrust point),
     vent_times.csv, lox_drain_times.csv, component_table.csv/.md, schematic.png and
     Hopper_Parameters.xlsx.

Monte Carlo (many runs in parallel): python Hopper_MonteCarlo.py

Run from this folder:  python Hopper.py
Plots can also be made later from the saved results without re-running:
  python Backend/Hopper_Plots.py list                      (everything that can be plotted)
  python Backend/Hopper_Plots.py plot ox_ladder engine     (named plots)
  python Backend/Hopper_Plots.py nodes "LOX Tank" COPV     (any nodes, elements, tracks...)
"""

import os
import time

import matplotlib

from HOPPER_DATA import *
from Backend.Hopper_Chamber_Map import HopperChamberMap
from Backend.Hopper_Thrust_Targeter import HopperThrustTargets
from Backend.Hopper_Prop_System import HopperPropSystem
from Backend.Hopper_Export import ExportEngineTargets, ExportActualPerformance




# Render plots off-screen unless they are to be shown (must precede any pyplot import).
if not PROP_SHOW_PLOTS:
    matplotlib.use("Agg")

from Backend.Hopper_Plots import HopperRun, make_plots




# ---- Chamber Map ---- #
if PROP_GENERATE_CHAMBER_MAP or not os.path.exists(CHAMBER_MAP_FILENAME):

    print("Generating chamber Pc/MR map...")

    HopperChamberMap(filename=CHAMBER_MAP_FILENAME)




# ---- Thrust Targets ---- #
# The flows OTV and FTV will be commanded to, one set per profile point.
print("Solving thrust targets...")

targets = HopperThrustTargets(
    PROP_THRUST_CURVE,
    chamber_map_filename=CHAMBER_MAP_FILENAME,
)

os.makedirs(PROP_RESULTS_DIRECTORY, exist_ok=True)
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

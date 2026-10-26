"""Run the Hopper vehicle fluids model (driver script, like Vespula.py).

Every input, including these run settings, is in HOPPER_DATA.py.

Steps, each switchable in HOPPER_DATA.py:
  1. Chamber map: build Hopper_Chamber_Map.h5 if asked to or if it is missing.
  2. Thrust targets: Pc and propellant flows for every thrust-profile point.
  3. Prop system: the FullFlow transient of the whole vehicle P&ID.
  4. Results: CSVs (vehicle_sim.py formats) and the plots listed in PROP_PLOTS
     in PROP_RESULTS_DIRECTORY.

Run from this folder:  python Hopper.py
Plots can also be made later from the saved results without re-running:
  python Hopper_Plots.py list                      (everything that can be plotted)
  python Hopper_Plots.py plot ox_ladder engine     (named plots)
  python Hopper_Plots.py nodes "LOX Tank" COPV     (any nodes, elements, tracks...)
"""

import os
import time

import matplotlib

from HOPPER_DATA import *
from Hopper_Chamber_Map import HopperChamberMap
from Hopper_Thrust_Targeter import HopperThrustTargets
from Hopper_Prop_System import HopperPropSystem
from Hopper_Export import ExportEngineTargets, ExportActualPerformance




# Render plots off-screen unless they are to be shown (must precede any pyplot import).
if not PROP_SHOW_PLOTS:
    matplotlib.use("Agg")

from Hopper_Plots import HopperRun, make_plots




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

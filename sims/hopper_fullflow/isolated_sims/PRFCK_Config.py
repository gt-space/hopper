"""PRFCK sizing inputs (Fuel tank pressurant check valve). Inputs only, no logic.

Section 1 is read from hopper_config.json through HOPPER_DATA: edit the JSON for those.
Sections 2 and 3 are this sizer's own inputs: edit them here.
Solver: Hopper_PRFCK_Sizing.py (same folder).
"""

import os
import sys
from pathlib import Path

import numpy as np

FULLFLOW_DIR = Path(__file__).resolve().parents[1]   # sims/hopper_fullflow
sys.path.insert(0, str(FULLFLOW_DIR))
os.chdir(FULLFLOW_DIR)

import HOPPER_DATA as D


# =====================================================================================
# 1. FROM hopper_config.json
# =====================================================================================

# ---- Run control ---- #
HISTORY_FILE = D.PROP_FILENAME            # full-flow transient history (Outputs/Hopper_Prop_System.h5)
TIME_STEP_S = D.PROP_DT                   # simulation.maximum_step_s, used if the full sim is run
SAVE_STEP_S = D.PROP_SAVE_DT              # simulation.save_interval_s
OUTPUT_DIR = Path(D.PROP_RESULTS_DIRECTORY) / "PRFCK_outputs"

# ---- This valve's leg and tank (flow_paths.fuel_pressurant, tanks.fuel) ---- #
VALVE = "PRFCK"
SIDE = "Fuel"
LEG = D.FUEL_PRESSURANT_LEG                           # GasLine, BangBangValve, Orifice, CheckValve
TANK = D.FUEL_TANK
TANK_RELIEF = "FRV"                # relief on this tank's ullage (side_paths)

# ---- Pressure sources, reliefs, fluids (bottles, side_paths, fluids) ---- #
PRESSURANT = D.PRESSURANT
AMBIENT_PA = D.AMBIENT_PRESSURE
COPV_INITIAL_PA = max(c.pressure for c in D.COPVS)
COPV_INITIAL_TEMPERATURE_K = max(c.temperature for c in D.COPVS)
COPV_MIN_TEMPERATURE_K = min(D.COPV_MIN_TEMPERATURES)
RELIEFS = {e.name: e for b in D.BRANCHES for e in b.elements if type(e).__name__ == "ReliefValve"}
COPV_RELIEF = RELIEFS["PRRV"]


# =====================================================================================
# 2. SIZING CRITERIA (this sizer's own)
# =====================================================================================

REUSE_HISTORY = True               # True: reuse HISTORY_FILE if it exists and matches the config;
                                   # False (or no/mismatched file): run the full flight sim (python Hopper.py)
REQUIRED_RATING_PSI = 5000.0       # methodology: pressure rating must be > this [psi]
RATING_FACTOR_ON_WORST_CASE = 1.0  # multiplies the static worst case (e.g. 1.5 for a proof factor)
CRACK_SWEEP_PSI = [1.0, 3.0, 10.0, 25.0]   # crack pressures to tabulate [psid]
CDA_LOSS_FRACTION = 0.05           # check may cut the leg effective CdA (bang-bang + orifice) by at most this
ORIFICE_GROWTH_FACTOR = 2.0        # also size for an orifice this much bigger than today's
CV_PER_IN2 = 38.0                  # Cv = 38 * CdA [in^2], the conversion used throughout hopper_config.json
CV_SWEEP = np.geomspace(0.01, 5.0, 300)    # Cv values for the min-Cv plot
FLOWING_THRESHOLD_KG_S = 1e-5      # leg mass flow above this counts as "check open and flowing" [kg/s]

# ---- Service requirements for this valve ---- #
OXYGEN_CLEAN_REQUIRED = False
SERVICE_REQUIREMENT = "RP-1 (hydrocarbon) compatible seals (outlet sees fuel ullage)"


# =====================================================================================
# 3. CANDIDATE VALVE (datasheet values of the part being evaluated; None = not entered yet)
# =====================================================================================

CANDIDATE_PART = None              # e.g. "Manufacturer model-number"
CANDIDATE_CDA_MM2 = 44.2           # flow area [mm^2]; None = use hopper_config.json's PRFCK CdA
CANDIDATE_CV = None                # Cv; if given, overrides CANDIDATE_CDA_MM2
CANDIDATE_RATED_PSI = None         # maximum working pressure rating [psi]
CANDIDATE_CRACK_PSI = None         # cracking pressure [psid]
CANDIDATE_MIN_TEMP_K = None        # lowest rated temperature [K]
CANDIDATE_MAX_TEMP_K = None        # highest rated temperature [K]
CANDIDATE_OXYGEN_CLEAN = None      # True / False: cleaned for oxygen service
CANDIDATE_SEAL_MATERIAL = None     # e.g. "PCTFE", "PTFE", "FKM (Viton)"
CANDIDATE_LEAK_RATING = None       # e.g. "bubble-tight at 600 psid"

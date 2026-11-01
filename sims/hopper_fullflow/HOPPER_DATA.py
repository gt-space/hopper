"""ALL THE INPUTS FOR THE HOPPER FULLFLOW MODEL ARE HERE.

Edit this file to change the vehicle, the thrust curve, the start-up sequence,
or how the model is run and plotted. No other file holds a vehicle number.
Then run python Hopper.py
"""

from pathlib import Path

from scipy.constants import *

from Backend.Hopper_Inputs import load_thrust_curve
from Backend.Hopper_Elements import (CLOSED, BangBangValve, Branch, CavitatingVenturi, CheckValve, Custom, GasLine, Injector,
                             Line, Orifice, Regulator, ReliefValve, Tank, ThrottleValve, Valve, tube_area, tube_id)




# =========================================================================== #
# ---- RUN SETTINGS ---- #
# =========================================================================== #
PROP_THRUST_CURVE_NAME     = "1dof"     # which entry of THRUST_CURVES (below) to fly
PROP_PROFILE_INTERPOLATION = "step"     # "step": hold each point until the next, as vehicle_sim.py does;
                                        # "linear": interpolate between points
PROP_DT                    = 0.05       # [s] time step. 0.05 s gives the same engine results as 0.01 s in
                                        # ~40 % of the time; the bang-bang is sampled less often, so tank
                                        # pressure overshoots ~10 psi more. Use 0.01 to study the tank band.
PROP_SAVE_DT               = 0.05       # [s] output interval in the HDF5 file
PROP_RTOL                  = 1e-6       # per-step residual tolerance
PROP_JACOBIAN_METHOD       = "2-point"  # "2-point" is ~2x faster than "3-point"
PROP_VERBOSE               = False      # True prints every solver iteration
PROP_STATISTICS            = True       # True prints each accepted step and a solver summary
PROP_SOLVE_PROP_SYSTEM     = True       # False re-exports results from an existing PROP_FILENAME
PROP_GENERATE_CHAMBER_MAP  = False      # True rebuilds the map (~2 min); it is also built if the file is missing

# Files
# Every path is anchored to this folder, so the scripts work from any working directory.
_HERE = Path(__file__).resolve().parent
CHAMBER_MAP_FILENAME = str(_HERE / "Backend" / "Hopper_Chamber_Map.h5")   # combustion-product map, built by Backend/Hopper_Chamber_Map.py
PROP_FILENAME = str(_HERE / "Outputs" / "Hopper_Prop_System.h5")          # FullFlow transient history

# Results
PROP_EXPORT_RESULTS        = True
PROP_RESULTS_DIRECTORY     = str(_HERE / "Outputs")
PROP_PLOTS                 = "all"      # "all", "deck" (the old sim's output deck), or a list of names,
                                        # e.g. ["nodes", "ox_ladder", "engine"]; see Hopper_Plots.PLOTS
PROP_PLOT_UNITS            = "E"        # "E" = psi / deg F (as the old sim), "SI" = Pa / K
PROP_PLOT_CSV              = False      # True also writes each plot's numbers to <plot>.csv
PROP_PLOT_PROFILE_TIME     = False      # True plots time from the start of the thrust profile (old sim's clock)
PROP_SHOW_PLOTS            = False      # True also opens the plot windows
PROP_REPORT                = True       # run summary, CdA table, vent times, component table, schematic, parameter workbook




# =========================================================================== #
# ---- THRUST CURVES ---- #
# =========================================================================== #
# Each entry is a CSV file (path relative to this folder) or an inline
# {profile time [s]: thrust [N]} dict. Select one with PROP_THRUST_CURVE_NAME.
# CSV format (e.g. a trajectory-sim export): a header row with a time column
# ("time [s]" or "[ms]") and a thrust column ("thrust [N]", "[kN]" or "[lbf]");
# other columns are ignored. Points closer than THRUST_CURVE_MIN_SPACING are
# resampled onto that spacing; use PROP_PROFILE_INTERPOLATION = "linear" for
# dense exports. See Hopper_Inputs.load_thrust_curve.
THRUST_CURVES = {
    "1dof": "thrust_curves/1dof.csv",   # 1dof trajectory sim profile (29.7 s), from Luca's vehicle_sim.py
}
THRUST_CURVE_MIN_SPACING = 0.1          # [s]
_curve = THRUST_CURVES[PROP_THRUST_CURVE_NAME]
PROP_THRUST_CURVE = (load_thrust_curve(_curve, THRUST_CURVE_MIN_SPACING) if isinstance(_curve, str) else _curve)
THRUST_CURVE_1DOF = load_thrust_curve(THRUST_CURVES["1dof"])   # used by Old_sims_comparison/compare_models.py




# =========================================================================== #
# ---- FLUIDS ---- #
# =========================================================================== #
OXIDIZER = "LOX"
FUEL = "RP-1"
PRESSURANT = "Nitrogen"  

AMBIENT_PRESSURE = 14.8 * psi




# =========================================================================== #
# ---- PRESSURANT: COPVs ---- #
# =========================================================================== #

COPV_COUNT = 3
COPV_UNIT_VOLUME = 9.0 * liter
COPV_VOLUME = COPV_COUNT * COPV_UNIT_VOLUME
PRESSURANT_INITIAL_PRESSURE = 4503.2 * psi   # [psi]
PRESSURANT_INITIAL_TEMPERATURE = 293.0       # [K]




# =========================================================================== #
# ---- TANKS ---- #
# =========================================================================== #
# volume [m^3], liquid_mass [kg], liquid_temperature [K]
# pressure [Pa] starting ullage pressure
# ullage_temperature [K]
# collapse_factor [-] COPV pressurant drawn / pressurant kept in the ullage (1 = no collapse)
LOX_TANK = Tank(
    "LOX Tank", propellant=OXIDIZER,
    volume=20.0 * liter, liquid_mass=19.43 , liquid_temperature=90.0,
    pressure=500.0 * psi, ullage_temperature=150.0,
    collapse_factor=2.0,   #EDIT assumed (user, 2026-10-27); Luca to supply conservative + Vespula values
    set_pressure=500.0 * psi, upper_limit=5.0 * psi, lower_limit=5.0 * psi,
)

FUEL_TANK = Tank(
    "Fuel Tank", propellant=FUEL,
    volume=16.0 * liter, liquid_mass=9.7, liquid_temperature=293.0,
    pressure=500.0 * psi, ullage_temperature=293.0,
    collapse_factor=1.1,   #EDIT assumed (user, 2026-10-27); Luca to supply conservative + Vespula values
    set_pressure=500.0 * psi, upper_limit=5.0 * psi, lower_limit=5.0 * psi,
)




# =========================================================================== #
# ---- PRESSURANT LEGS (COPV -> ... -> tank ullage) ---- #
# =========================================================================== #
TUBE_ROUGHNESS = 1.5e-5   # [m] vehicle_sim.py Line roughness, all tubes

LOX_PRESSURANT_LEG = [
    GasLine("LOX Pressurant Line", inner_diameter=4.92 * milli, length=1.0, #meter, EDIT
            roughness=TUBE_ROUGHNESS, outlet_node="OBB Inlet"),
    BangBangValve("OBB", cda=18.68 * milli**2),
    Orifice("OR-1", cda=2 * 1.02 * milli**2),              #EDIT
    CheckValve("PROCK", cda=12.0 * milli**2),              #EDIT
]

FUEL_PRESSURANT_LEG = [
    GasLine("Fuel Pressurant Line", inner_diameter=4.92 * milli, length=1.0, #meter, EDIT
            roughness=TUBE_ROUGHNESS, outlet_node="FBB Inlet"),
    BangBangValve("FBB", cda=18.68 * milli**2),             
    Orifice("OR-2", cda=1.02 * milli**2),                   #EDIT
    CheckValve("PRFCK", cda=12.0 * milli**2),               #EDIT
]




# =========================================================================== #
# ---- FEED LINES (tank liquid -> ... -> main chamber) ---- #
# =========================================================================== #
LOX_TUBE_ID = tube_id(3 / 8 * inch, 0.028 * inch)    # EDIT mm
FUEL_TUBE_ID = tube_id(1 / 4 * inch, 0.028 * inch)   # EDIT mm

INJECTOR_PASSAGE_LENGTH = 0.01   # EDIT [m]
REGEN_PASSAGE_LENGTH = 0.30      # EDIT [m]

LOX_FEED = [
    Line("LOX Tank to OTV Line", inner_diameter=LOX_TUBE_ID, length=1.0, roughness=TUBE_ROUGHNESS,
         outlet_node="OTV Inlet"),
    ThrottleValve("OTV", max_cda=0.9 * tube_area(LOX_TUBE_ID),          #EDIT
                  full_scale_mass_flow=0.777, outlet_node="OTV Outlet"),  
    Line("OTV to LOX Injector Line", inner_diameter=LOX_TUBE_ID, length=0.2032, roughness=TUBE_ROUGHNESS,
         outlet_node="LOX Injector Manifold"),
    Injector("OINJ", cda=15.3 * milli**2, length=INJECTOR_PASSAGE_LENGTH),
]

FUEL_FEED = [
    Line("Fuel Tank to FTV Line", inner_diameter=FUEL_TUBE_ID, length=0.3048, roughness=TUBE_ROUGHNESS,
         outlet_node="FTV Inlet"),
    ThrottleValve("FTV", max_cda=0.9 * tube_area(FUEL_TUBE_ID),         # EDIT: 90 % of tube area (17.2 mm^2)
                  full_scale_mass_flow=0.388, outlet_node="FTV Outlet"),
    Line("FTV to Regen Line", inner_diameter=FUEL_TUBE_ID, length=0.2032, roughness=TUBE_ROUGHNESS,
         outlet_node="Regen Inlet Manifold"),
    Orifice("Regen Jacket", cda=25.0 * milli**2, length=REGEN_PASSAGE_LENGTH,   # EDIT
            outlet_node="Fuel Injector Manifold"),
    Injector("FINJ", cda=10.6 * milli**2, length=INJECTOR_PASSAGE_LENGTH),
]

# Example additions (uncomment to use):
#   Valve("OMV", cda=40 * milli**2, schedule=([0.0, 0.2, 0.3], [0.0, 0.0, 1.0]))   # scheduled main valve
#   CavitatingVenturi("FCV", throat_area=6 * milli**2)                               # FullFlow cavitating venturi
#   Custom("My Part", build_function)                                                 # any FullFlow component




# =========================================================================== #
# ---- VENTS, RELIEFS, PURGES, DRAINS (P&ID side paths) ---- #
# =========================================================================== #
# Branch(name, source, [elements in series], sink)
#   source: "COPV", "LOX Ullage", "Fuel Ullage", "LOX Tank Liquid", "Fuel Tank Liquid",
#           a feed junction ("OTV Inlet"), or a regulator outlet ("PR-REG Outlet")
#   sink:   "ambient" (default), "Main Chamber", or a feed junction ("OTV Outlet")
#   into_sink=True: a liquid branch whose flow goes INTO the main chamber (the torch)
# Valve schedules are (times [s], openings 0-1); CLOSED keeps a valve shut all run.
# ReliefValve crack / full-open pressures are gauge [Pa] (full open defaults to 1.1 x crack).
# sheet = Hopper Fluid Components Mastersheet, info = Extra Component Info
REGULATORS = [
    Regulator("PR-REG", source="COPV", set_pressure=500.0 * psi, cda=5.0 * milli**2),   #EDIT GUESS: not in sheet
]

# Torch igniter: OIMV taps the LOX feed above OTV, FIMV the fuel feed above FTV
# (P&ID engine page); both feed the torch, which fires into the main chamber.
TORCH_OPEN_TIME = 0.20         #EDIT GUESS [s] OIMV + FIMV open (main valves start at MAIN_VALVE_OPEN_TIME)
TORCH_DURATION = 0.20          #EDIT [s] torch on-time (~0.2 s)
TORCH_VALVE_STROKE = 0.01      #EDIT GUESS [s] OIMV / FIMV opening and closing time
TORCH_SCHEDULE = ([0.0, TORCH_OPEN_TIME, TORCH_OPEN_TIME + TORCH_VALVE_STROKE,
                   TORCH_OPEN_TIME + TORCH_DURATION - TORCH_VALVE_STROKE, TORCH_OPEN_TIME + TORCH_DURATION],
                  [0.0, 0.0, 1.0, 1.0, 0.0])

JAKSA_CRYO_SOLENOID_CDA = 9.4 * milli**2   # Jaksa 2/2 cryo solenoid (info), CdA from sheet OVNT tab; used for OVNT and OFD

BRANCHES = [
    # Pressurant manifold
    Branch("PRVNT Vent", "COPV", [
        Valve("PRVNT", cda=1.07 * milli**2, schedule=CLOSED)]),       # Marotta PLV34 (info); 1x PLV34 = 1.07 mm^2 (sheet FRV tab)
    Branch("PRRV Relief", "COPV", [
        ReliefValve("PRRV", cda=5.0 * milli**2, crack_pressure=4950.0 * psi)]),   #EDIT GUESS: Swagelok SS-4R3A (info), sheet tab empty
    Branch("ORV-1 Relief", "PR-REG Outlet", [
        ReliefValve("ORV-1", cda=2.0 * milli**2, crack_pressure=600.0 * psi)]),   #EDIT GUESS
    # GN2 engine purges, fed from PR-REG
    Branch("LOX Purge", "PR-REG Outlet", [
        Valve("OPV", cda=2.0 * milli**2, schedule=CLOSED),            #EDIT GUESS CdA
        CheckValve("CV-5", cda=5.0 * milli**2)],                      #EDIT GUESS CdA
        sink="OTV Outlet"),                                           # P&ID: joins the LOX line below OTV
    Branch("Fuel Purge", "PR-REG Outlet", [
        Valve("FPV", cda=2.0 * milli**2, schedule=CLOSED),            #EDIT GUESS CdA
        CheckValve("CV-4", cda=5.0 * milli**2)],                      #EDIT GUESS CdA
        sink="FTV Outlet"),                                           # P&ID: joins the fuel line below FTV
    Branch("Torch Fuel Purge", "PR-REG Outlet", [
        Valve("TFPV", cda=2.0 * milli**2, schedule=CLOSED),           #EDIT GUESS CdA
        CheckValve("CV-3", cda=5.0 * milli**2)],                      #EDIT GUESS CdA
        sink="Main Chamber"),                                         # P&ID: torch fuel line below FIMV -> torch
    # LOX tank
    Branch("OVNT Vent", "LOX Ullage", [
        Valve("OVNT", cda=JAKSA_CRYO_SOLENOID_CDA, schedule=CLOSED)]), # Jaksa 2/2 cryo solenoid (info), 9.4 mm^2 (sheet OVNT)
    Branch("ORV Relief", "LOX Ullage", [
        ReliefValve("ORV", cda=30.8 * milli**2, crack_pressure=550.0 * psi)]),    # RegO PRV9432T550 (info); target CdA (sheet ORV)
    Branch("OFD Drain", "LOX Tank Liquid", [
        Valve("OFD", cda=JAKSA_CRYO_SOLENOID_CDA, schedule=CLOSED)]), # same valve as OVNT (Jaksa 2/2 cryo solenoid)
    # Fuel tank
    Branch("FVNT Vent", "Fuel Ullage", [
        Valve("FVNT", cda=15.0 * milli**2, schedule=CLOSED)]),        # UNOX 3/8 NO solenoid (info); target CdA (sheet FVNT)
    Branch("FRV Relief", "Fuel Ullage", [
        ReliefValve("FRV", cda=14.67 * milli**2, crack_pressure=550.0 * psi)]),   #EDIT PRV9432T550 = 550 psig (info); sheet FRV says 525
]

# Torch igniter branches (flow is delivered into the main chamber). Off by default:
# it fires ~0.2 s and barely changes the burn; True adds it to BRANCHES.
TORCH_ENABLED = False
TORCH_BRANCHES = [
    Branch("Torch Ox", "OTV Inlet", [
        Valve("OIMV", cda=2.0 * milli**2, schedule=TORCH_SCHEDULE),   #EDIT GUESS CdA
        Orifice("HX-1", cda=2.0 * milli**2),                          #EDIT GUESS CdA
        Orifice("Torch Ox Injector", cda=0.08 * milli**2)],           #EDIT GUESS CdA (~7 g/s LOX at 500 psi)
        sink="Main Chamber", into_sink=True),
    Branch("Torch Fuel", "FTV Inlet", [
        Valve("FIMV", cda=2.0 * milli**2, schedule=TORCH_SCHEDULE),   #EDIT GUESS CdA
        Orifice("Torch Fuel Injector", cda=0.06 * milli**2)],         #EDIT GUESS CdA (~4.5 g/s RP-1, torch MR ~1.6)
        sink="Main Chamber", into_sink=True),
]
if TORCH_ENABLED:
    BRANCHES += TORCH_BRANCHES

# Vent-time analysis (Hopper_Vent.py): each commanded vent off a tank ullage
# blows the tank down from set + upper limit to ambient + VENT_END_GAUGE_PRESSURE.
VENT_END_GAUGE_PRESSURE = 1.0 * psi   #EDIT GUESS: "vented" = 1 psig
VENT_TIME_REQUIREMENTS = {            # [s] by branch name
    "OVNT Vent": 60.0,                # sheet OVNT target vent time
    "FVNT Vent": 30.0,                # sheet FVNT vent time
}
VENT_DT = 0.1                         # [s]
VENT_MAX_TIME = 300.0                 # [s]

# LOX drain calculator (Hopper_Drain.py): time to empty the LOX tank through the
# OFD Drain branch (+ LOX_DRAIN_EXTRA_ELEMENTS) with the ullage held at each gauge
# pressure in LOX_DRAIN_ULLAGE_PRESSURES (0 = vented, gravity only).
LOX_TANK_HEIGHT = 0.6167                  # [m] sheet ORV tab; cross-section = LOX_TANK.volume / height
LOX_DRAIN_BRANCH = "OFD Drain"
LOX_DRAIN_EXTRA_ELEMENTS = [              #EDIT ground side of the drain (OQD, GSE line) once known, e.g.
    # Orifice("OQD", cda=20.0 * milli**2),
]
LOX_DRAIN_OUTLET_DROP = 0.0               #EDIT GUESS [m] drain exit below the tank bottom
LOX_DRAIN_ULLAGE_PRESSURES = [0.0, 10.0 * psi, 25.0 * psi, 50.0 * psi]   #EDIT GUESS [Pa gauge]
LOX_DRAIN_END_FRACTION = 0.01             # "drained" = 1 % of the load left
LOX_DRAIN_TIME_REQUIREMENT = None         #EDIT [s], None = no requirement
DRAIN_STEPS = 400                         # time steps per drain run (length set from an estimate of the drain time)
DRAIN_MAX_TIME = 3600.0                   # [s]

# P&ID hardware with no flow path in a burn: in the component table and schematic only.
NOT_MODELLED = {
    "PRQD": "Pressurant fill quick disconnect (disconnected in flight)",
    "OQD": "LOX fill quick disconnect below OFD (disconnected in flight)",
    "FQD": "Fuel fill quick disconnect, McMaster 52495K63 (disconnected in flight)",
    "ODiff": "Oxygen diffuser on the OVNT exhaust (sizing)",
    "Torch spark igniter": "Spark / exciter; torch combustion itself is not modelled (its propellant goes into the main chamber)",
    "PRPT-1 and sensors": "Instrumentation",
}
if not TORCH_ENABLED:
    NOT_MODELLED["OIMV / FIMV / HX-1"] = "Torch igniter feeds, switched off (TORCH_ENABLED = False)"




# =========================================================================== #
# ---- ENGINE ---- #
# =========================================================================== #
MIXTURE_RATIO = 1.8               # target O/F
THROAT_AREA = 0.00067226          # [m^2]
EXIT_AREA = 0.00235559904         # [m^2]
EXPANSION_RATIO = EXIT_AREA / THROAT_AREA   # 3.504

# Efficiency vs chamber pressure. Same parameters as Luca's previous sim; linear
# ramp from LOW at PC_LOW to HIGH at PC_LOW + PC_SPAN, held flat outside it:
#   eta = LOW + (HIGH - LOW) * clip((Pc_psia - PC_LOW) / PC_SPAN, 0, 1)
# (eta_c* 0.80 at 100 psia -> 0.85 at 300 psia; eta_cf 0.94 at 100 -> 0.97 at 250.)
CSTAR_EFFICIENCY_LOW = 0.80
CSTAR_EFFICIENCY_HIGH = 0.85
CSTAR_EFFICIENCY_PC_LOW = 100.0     # [psia]
CSTAR_EFFICIENCY_PC_SPAN = 200.0    # [psi]
CF_EFFICIENCY_LOW = 0.94
CF_EFFICIENCY_HIGH = 0.97
CF_EFFICIENCY_PC_LOW = 100.0        # [psia]
CF_EFFICIENCY_PC_SPAN = 150.0       # [psi]


NOMINAL_CHAMBER_PRESSURE = 300.0 * psi

# Combustion map grid (rebuild with PROP_GENERATE_CHAMBER_MAP = True after a
# change). The low end reaches below ambient so start-up and tail-off stay
# inside the map.
CHAMBER_MAP_PRESSURE_LOWER_LIMIT = 10.0 * psi
CHAMBER_MAP_PRESSURE_UPPER_LIMIT = 610.0 * psi
CHAMBER_MAP_PRESSURE_POINTS = 31
CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT = 0.5
CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT = 4.0
CHAMBER_MAP_MIXTURE_RATIO_POINTS = 36




# =========================================================================== #
# ---- SEQUENCE ---- #
# =========================================================================== #
# GUESS (vehicle_sim.py started at full flow at t = 0): OTV/FTV flow targets
# ramp up, the thrust profile starts when they finish (0.5 s), they ramp to
# zero after the last profile point, then the tail-off is simulated. The torch
# igniter is not modelled; propellant reaching the chamber ignites at once.
MAIN_VALVE_OPEN_TIME = 0.30   # [s] OTV + FTV start opening
MAIN_VALVE_RAMP_TIME = 0.20   # [s] OTV/FTV opening (and closing) time
THROTTLE_STEP_TIME = 0.01     # [s] OTV/FTV travel time for one "step" change between profile points
RAMP_BREAKPOINT_SPACING = 0.01   # [s] solver steps no longer than this through the opening / closing ramps
TAIL_OFF_TIME = 0.30          # [s] simulated after OTV/FTV have closed




# =========================================================================== #
# ---- MONTE CARLO (Hopper_MonteCarlo.py) ---- #
# =========================================================================== #
MC_RUNS = 100                 # run 0 is nominal
MC_SEED = 1                   # same seed = same draws
MC_WORKERS = 16               # parallel processes (this PC has 20 cores; each run is ~1 min)
MC_SAVE_DT = 0.25             # [s] output interval for each run's file
MC_KEEP_FILES = False         # True keeps every run's HDF5 file in MC_DIRECTORY
MC_DIRECTORY = str(_HERE / "Outputs" / "MonteCarlo")
MC_THRUST_CONTROL = "thrust"  # "thrust": controller holds the thrust target with each run's engine;
                              # "flow": every run gets the nominal flow commands
# {group: (distribution, size, "rel" or "abs", inputs, "each" or "shared")}
#   distribution  "uniform": size = half-width (every value in +-size equally likely);
#                 "normal": size = 1 sigma, clipped at 3 sigma
#   rel / abs     size is a fraction of each input's nominal / in the input's own (SI) units
#   inputs        names above, "NAME.attribute" (e.g. "LOX_TANK.pressure"), "<element>.attribute"
#                 (e.g. "OINJ.cda_value"), or MC_ALL_CDAS = every CdA that carries flow in a run
#                 (pressurant legs, feeds incl. OTV/FTV max CdA, injectors, regen, reliefs,
#                 regulator; never-opened vents / purges / drains are left out)
#   each / shared "each": every input gets its own draw; "shared": one draw moves them all
#                 together (e.g. the whole c* efficiency curve up or down)
MC_ALL_CDAS = "ALL_CDAS"
MC_DISPERSIONS = {
    "CdA": ("uniform", 0.10, "rel", MC_ALL_CDAS, "each"),                        #EDIT +-10 % (user)
    "Initial tank pressure": ("uniform", 0.02, "rel",
                              ["LOX_TANK.pressure", "FUEL_TANK.pressure"], "each"),  #EDIT +-2 % (user)
    "Initial COPV pressure": ("uniform", 0.02, "rel", ["PRESSURANT_INITIAL_PRESSURE"], "each"),  #EDIT +-2 % (user)
    "Propellant load": ("uniform", 0.01, "rel",
                        ["LOX_TANK.liquid_mass", "FUEL_TANK.liquid_mass"], "each"),  #EDIT +-1 % (user)
    "Collapse factor": ("uniform", 0.3, "abs",
                        ["LOX_TANK.collapse_factor", "FUEL_TANK.collapse_factor"], "each"),  #EDIT +-0.3 (user)
    "c* efficiency": ("uniform", 0.02, "rel",
                      ["CSTAR_EFFICIENCY_LOW", "CSTAR_EFFICIENCY_HIGH"], "shared"),  #EDIT +-2 % (user)
    "Cf efficiency": ("uniform", 0.01, "rel",
                      ["CF_EFFICIENCY_LOW", "CF_EFFICIENCY_HIGH"], "shared"),    #EDIT +-1 % (user)
    "COPV temperature": ("uniform", 10.0, "abs", ["PRESSURANT_INITIAL_TEMPERATURE"], "each"),  #EDIT GUESS +-10 K
    "Mixture ratio": ("uniform", 0.03, "rel", ["MIXTURE_RATIO"], "each"),        #EDIT GUESS +-3 % controller MR bias
}
# Physical limits applied after a draw, by attribute name: (min, max), None = no limit.
MC_INPUT_LIMITS = {
    "collapse_factor": (1.0, None),   # below 1 the ullage would gain gas from nowhere
}




# =========================================================================== #
# ---- PLOTS ---- #
# =========================================================================== #
# Flow-separation line on the exit-pressure plot: exit pressure below this
# multiple of ambient risks separation (the old deck's "0.4 x Pamb" line).
NOZZLE_SEPARATION_RATIO = 0.4

# Run report (Hopper_Report.py)
REPORT_THRUST_ON_FRACTION = 0.10   # burn time = time thrust is above this fraction of the max target
REPORT_THRUST_TOLERANCE = 0.02     # each profile point's thrust must be within this fraction of target




# =========================================================================== #
# ---- SOLVER NUMERICS (rarely need changing) ---- #
# =========================================================================== #
# Mixture ratio is (ox + MR x floor) / (fuel + floor), so it stays defined
# while nothing flows; negligible once propellant moves.
MIXTURE_RATIO_FLOW_FLOOR = 1.0e-4   # [kg/s]
# Within this much of ambient, nozzle flow is made linear in (Pc - ambient),
# so the solver has a usable slope at zero flow (see Hopper_Components.py).
NOZZLE_LINEAR_PRESSURE_DROP = 0.5 * psi

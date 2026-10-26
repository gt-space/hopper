"""ALL THE INPUTS FOR THE HOPPER FULLFLOW MODEL ARE HERE.

Edit this file to change the vehicle, the thrust curve, the start-up sequence,
or how the model is run and plotted. No other file holds a vehicle number.
Then run python Hopper.py
"""

from scipy.constants import *

from Hopper_Elements import (BangBangValve, CavitatingVenturi, CheckValve, Custom, GasLine, Injector, Line,
                             Orifice, Tank, ThrottleValve, Valve, tube_area, tube_id)




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
CHAMBER_MAP_FILENAME = "Hopper_Chamber_Map.h5"   # combustion-product map, built by Hopper_Chamber_Map.py
PROP_FILENAME = "Hopper_Prop_System.h5"          # FullFlow transient history

# Results
PROP_EXPORT_RESULTS        = True
PROP_RESULTS_DIRECTORY     = "Hopper_Results"
PROP_PLOTS                 = "all"      # "all", "deck" (the old sim's output deck), or a list of names,
                                        # e.g. ["nodes", "ox_ladder", "engine"]; see Hopper_Plots.PLOTS
PROP_PLOT_UNITS            = "E"        # "E" = psi / deg F (as the old sim), "SI" = Pa / K
PROP_PLOT_CSV              = False      # True also writes each plot's numbers to <plot>.csv
PROP_PLOT_PROFILE_TIME     = False      # True plots time from the start of the thrust profile (old sim's clock)
PROP_SHOW_PLOTS            = False      # True also opens the plot windows




# =========================================================================== #
# ---- THRUST CURVES ---- #
# =========================================================================== #
# Each curve maps time from the start of the profile [s] to thrust [N]. Add a
# curve here and select it with PROP_THRUST_CURVE_NAME above. Points must be
# at least 0.1 s apart.
THRUST_CURVES = {
    # The full flight profile from the 1dof trajectory sim (29.7 s), from vehicle_sim.py.
    "1dof": {
        0.0: 1115.9, 0.1: 1176.9, 0.3: 1238.0, 0.6: 1309.1, 1.0: 1365.7,
        1.6: 1382.5, 2.5: 1347.9, 3.5: 1276.7, 4.5: 1209.9, 5.5: 1007.4,
        6.5: 890.0,  7.5: 890.0,  8.5: 890.0,  9.5: 890.0,  10.5: 970.4,
        11.5: 1015.4, 12.5: 1034.7, 13.5: 1040.9, 14.5: 1040.6, 15.5: 896.9,
        16.3: 890.0, 17.3: 890.0, 18.3: 890.0, 19.3: 890.0, 20.3: 890.0,
        21.3: 890.0, 21.8: 1074.1, 22.4: 1213.5, 23.1: 1218.3, 24.1: 1172.4,
        25.1: 1115.0, 26.1: 1066.2, 27.1: 1030.0, 28.1: 1004.8, 29.1: 987.6, 29.7: 979.9,
    },
}
THRUST_CURVE_1DOF = THRUST_CURVES["1dof"]   # kept for comparison/compare_models.py
PROP_THRUST_CURVE = THRUST_CURVES[PROP_THRUST_CURVE_NAME]




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
# radius [m] sets the liquid free-surface area for ullage heat transfer
# ullage_htc [W/m^2/K] 
LOX_TANK = Tank(
    "LOX Tank", propellant=OXIDIZER,
    volume=20.0 * liter, liquid_mass=19.43 , liquid_temperature=90.0,
    pressure=500.0 * psi, ullage_temperature=150.0,
    radius=0.2, ullage_htc=500.0,
    set_pressure=500.0 * psi, upper_limit=5.0 * psi, lower_limit=5.0 * psi,
)

FUEL_TANK = Tank(
    "Fuel Tank", propellant=FUEL,
    volume=16.0 * liter, liquid_mass=9.7, liquid_temperature=293.0,
    pressure=500.0 * psi, ullage_temperature=293.0,
    radius=0.1, ullage_htc=0.0,
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
# ---- ENGINE ---- #
# =========================================================================== #
MIXTURE_RATIO = 1.8               # target O/F
THROAT_AREA = 0.00067226          # [m^2]
EXIT_AREA = 0.00235559904         # [m^2]
EXPANSION_RATIO = EXIT_AREA / THROAT_AREA   # 3.504

# Efficiency vs chamber pressure. Same parameters as Luca's previous sim.
#   eta = LOW + (HIGH - LOW) * max(1, (Pc_psia - PC_LOW) / PC_SPAN)
# (max(1, ...) holds eta_c* at 0.85 up to 300 psia and eta_cf at 0.97 up to
# 250 psia; above those they keep rising with Pc.)
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
TAIL_OFF_TIME = 0.30          # [s] simulated after OTV/FTV have closed




# =========================================================================== #
# ---- PLOTS ---- #
# =========================================================================== #
# Flow-separation line on the exit-pressure plot: exit pressure below this
# multiple of ambient risks separation (the old deck's "0.4 x Pamb" line).
NOZZLE_SEPARATION_RATIO = 0.4




# =========================================================================== #
# ---- SOLVER NUMERICS (rarely need changing) ---- #
# =========================================================================== #
# Mixture ratio is (ox + MR x floor) / (fuel + floor), so it stays defined
# while nothing flows; negligible once propellant moves.
MIXTURE_RATIO_FLOW_FLOOR = 1.0e-4   # [kg/s]
# Within this much of ambient, nozzle flow is made linear in (Pc - ambient),
# so the solver has a usable slope at zero flow (see Hopper_Components.py).
NOZZLE_LINEAR_PRESSURE_DROP = 0.5 * psi

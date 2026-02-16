import time
import math
import numpy as np
import general_fluid_network as gfn
from general_fluid_network import Node, Ambient, Tank, BangBang, ThrottleValve, Line, Network

# ==============================================================================
# 1. CONSTANTS & GEOMETRY
# ==============================================================================
PSI_TO_PA = 6894.75729
IN_TO_M = 0.0254
MM2_TO_M2 = 1e-6

# --- Tube Geometry Calculations ---
# Fuel Side: AL6061 (OD 0.375", Thk 0.028")
fuel_od = 0.375 * IN_TO_M
fuel_thk = 0.028 * IN_TO_M
fuel_id = fuel_od - (2 * fuel_thk)

# Ox Side: Stainless 304 (OD 0.375", Thk 0.01")
ox_od = 0.375 * IN_TO_M
ox_thk = 0.01 * IN_TO_M
ox_id = ox_od - (2 * ox_thk)

# Flow Restrictions
cda_orifice = 1.02 * MM2_TO_M2  # Restriction for regulation
# CdA for main valves (Approximate based on line size or user intent)
cda_valve_max = np.pi * (ox_id/2)**2 

# ==============================================================================
# 2. NODES (TANKS & AMBIENT)
# ==============================================================================

# Pressurant Supply
copv = Node("Nitrogen", m=4.8, V=15.47, T=293, name="COPV")

# Oxidizer Tank (LOX)
ox_tank = Tank(V_total_L=20, 
               fluid_liq='Oxygen', 
               m_liq=19.43, 
               T_liq=90, 
               fluid_ullage="Nitrogen", 
               P_ullage=500*PSI_TO_PA, 
               T_ullage=150,
               radius=0.2, 
               htc=500,
               name="OxTank")

# Fuel Tank (n-Dodecane)
fu_tank = Tank(V_total_L=16, 
               fluid_liq='n-Dodecane', 
               m_liq=9.7, 
               T_liq=293, 
               fluid_ullage="Nitrogen", 
               P_ullage=500*PSI_TO_PA, 
               T_ullage=293,
               radius=0.1, 
               htc=0,
               name="FuTank")

# Ambient / Overboard / Engine
amb = Ambient(fluid="Air", P=14.7*PSI_TO_PA, T=293.15, name="Ambient")

# ==============================================================================
# 3. PRESSURIZATION SYSTEM (Bang-Bang Regulators)
# ==============================================================================

# COPV -> Tank Ullages
# Note: Connects directly to the Tank object; the network automatically routes to ullage
obb = BangBang(CdA=cda_orifice, 
               target_pressure=(500*PSI_TO_PA),
               hysteresis=(5*PSI_TO_PA),
               name="OxReg")

fbb = BangBang(CdA=cda_orifice, 
               target_pressure=(500*PSI_TO_PA),
               hysteresis=(5*PSI_TO_PA),
               name="FuReg")

# ==============================================================================
# 4. FEED SYSTEMS (Lines & Valves)
# ==============================================================================

# --- Oxidizer Feed String ---
# Path: Tank -> Line1 -> Valve -> Line2 -> Ambient
ox_line1 = Line(ID=ox_id, length=0.7735, name="OxLine_Tank2OTV")
# Setting normal_state here sets the initial condition at T=0
otv      = ThrottleValve(CdA_max=cda_valve_max, name="OxThrottle", normal_state=0.36) 
ox_line2 = Line(ID=ox_id, length=0.1524, name="OxLine_OTV2TCA")

# Stitch them together (creates intermediate Junction nodes automatically)
ox_conns, ox_nodes = gfn.connect_series(ox_tank, amb, [ox_line1, otv, ox_line2])

# --- Fuel Feed String ---
# Path: Tank -> Line1 -> Valve -> Line2 -> Ambient
fu_line1 = Line(ID=fuel_id, length=0.3048, name="FuLine_Tank2FTV")
ftv      = ThrottleValve(CdA_max=cda_valve_max, name="FuThrottle", normal_state=0.19)
fu_line2 = Line(ID=fuel_id, length=0.3302, name="FuLine_FTV2TCA")

# Stitch them together
fu_conns, fu_nodes = gfn.connect_series(fu_tank, amb, [fu_line1, ftv, fu_line2])

# ==============================================================================
# 5. NETWORK ASSEMBLY
# ==============================================================================

# Start with the Pressurization graph
graph = {
    obb: (copv, ox_tank),
    fbb: (copv, fu_tank)
}

# Add the Oxidizer Feed String
for comp, n1, n2 in ox_conns:
    graph[comp] = (n1, n2)

# Add the Fuel Feed String
for comp, n1, n2 in fu_conns:
    graph[comp] = (n1, n2)

network = Network(graph)

# ==============================================================================
# 6. SIMULATION CONTROL
# ==============================================================================

# Throttle Profile (Time -> Valve State 0-1)
# Note: keys must handle float precision (sim uses round(t, 4))
actions = {
    # 0.0 state is handled by 'normal_state' in init above
    0.1: (otv, 0.38),
    0.2: (ftv, 0.19),
    0.3: (otv, 0.4),
    0.4: (ftv, 0.2),
    0.6: (otv, 0.42),
    0.7: (ftv, 0.21),
    1.0: (otv, 0.44),
    1.1: (ftv, 0.22),
    1.6: (otv, 0.45),
    1.7: (ftv, 0.22),
    2.5: (otv, 0.44),
    2.6: (ftv, 0.22),
    3.5: (otv, 0.41),
    3.6: (ftv, 0.21),
    4.5: (otv, 0.39),
    4.6: (ftv, 0.2),
    5.5: (otv, 0.33),
    5.6: (ftv, 0.16),
    6.5: (otv, 0.29),
    6.6: (ftv, 0.14),
    # ... fill in rest ...
    29.7: (otv, 0.32),
    29.8: (ftv, 0.16)
}

print("Starting Simulation...")
print(f"Initial State: COPV={copv.P/PSI_TO_PA:.0f} psi, OxTank={ox_tank.P/PSI_TO_PA:.0f} psi")

dt = 0.1  # Timestep
runtime = 25.0 # Set to desired length
start_time = time.time()

network.sim(runtime, dt, actions, verbose_steps=100)

print(f"Simulation finished in {time.time() - start_time:.2f}s")

# ==============================================================================
# 7. PLOTTING
# ==============================================================================

# Nodes of interest (Includes the new auto-generated junctions!)
plot_nodes = [copv, ox_tank, fu_tank] + ox_nodes + fu_nodes

# Connections of interest
plot_conns = [obb, fbb, ox_line1, otv, ox_line2, fu_line1, ftv, fu_line2]

network.plot_nodes_overlay(plot_nodes, title="System Pressures & Masses", units="E")
network.plot_connections_overlay(plot_conns, title="Flow Rates & Line Dynamics", units="E")
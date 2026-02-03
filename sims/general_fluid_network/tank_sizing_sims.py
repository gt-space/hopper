from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network, Regulator, Tank, BangBang
#################### TEST CONFIGS #########################
import math


# GAS PULL
# test_node2 = Node("N2O", 35, 60, 293, "Gas Pull")
# amb_node = Ambient()
# test_connection1 = Connection(0.0000063, 0, 0)
# test_connection2 = Connection(0.00001521, 0, 1)
# test_connection3 = Connection(0.0000001742, 0, 0)
# test_network = Network({test_connection1: (test_node2, amb_node)})
# test_network.sim(1000, 1)
# test_network.plot_nodes_overlay((test_node1, test_node2), units="E")

# FILL + FULL THROTTLE
# target mass: 40 kg
# n = 7 # num bottles
# vol = 50 # L
# bottle_temp = 305 # K
# burn_duration = 55 # s
# fill_duration = 1100 # s
# vehicle_tank = Node("N2O", 1, vol, 293, "Vehicle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, bottle_temp, "Fill Tanks")
# amb_node = Ambient()
# chamber = Ambient(P=101325*350*1.25/14.7)
# fill_line = Connection(0.00005, 0, 0)
# vent_line = Connection(0.0000005, 0, 1, False)
# tv1 = ThrottleValve(1, target_mdot=0.0, normal_state=0.0)
# test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), tv1: (vehicle_tank, chamber)})
# test_network.sim(fill_duration+burn_duration, 1, {150: (vent_line, True), fill_duration-10: (fill_line, False), fill_duration-1: (vent_line, False), fill_duration: (tv1, 0.57*0.5), fill_duration+10: (tv1, 0.57)})
# test_network.plot_nodes_overlay((fill_tanks, vehicle_tank), title=f"{burn_duration}s, {vol}L, {n} Bottles, {bottle_temp}K", units="E")
# test_network.plot_connections_overlay([tv1], units="E")

# TV SIM
# vehicle_tank1 = Node("N2O", 40, 50, 288, "1")
# vehicle_tank2 = Node("N2O", 40, 50, 288, "2")

# chamber = Ambient(P=101325*435/14.7)
# tv1 = ThrottleValve(1, target_mdot=0.57, normal_state=0.57, name="1")
# tv2 = ThrottleValve(1, target_mdot=0.57, normal_state=0.57, name="2")

# test_network = Network({tv1: (vehicle_tank1, chamber), tv2: (vehicle_tank2, chamber)})
# throttle_profile = {30: (tv2, 0.57/2)}
# test_network.sim(57, 0.1, throttle_profile)
# test_network.plot_nodes_overlay([vehicle_tank1, vehicle_tank2], units="E")
# test_network.plot_connections_overlay([tv1,tv2], units="E")

# LIQUID PULL SMALL TANK
# n = 1
# vehicle_tank = Node("N2O", 36, 50, 293, "Vehicle Tank")
# vt2 = Node("N2O", 9.0718 * n, 13.4 * n, 293, "Bottle")
# amb_node = Ambient(P=101325*300/14.7)
# fluid_system = Connection(0.000009, 0, 0)
# fluid_system2 = Connection(0.000009, 0, 0)
# test_network = Network({fluid_system: (vehicle_tank, amb_node), fluid_system2: (vt2, amb_node)})
# test_network.sim(80, 1)
# test_network.plot_nodes_overlay((vehicle_tank, vt2), units="E")


# SUBCOOLED COPV 
# copv = Node("Nitrogen", 2.4, 6.61, 293, "COPV")
# vehicle_tank = Node("N2O", 36, 60, 287, "Liquid Pull")
# amb_node = Ambient()
# tank_substitute = Ambient(fluid="Nitrogen", P=101325*400/14.7)
# reg = Regulator(0.00000063, 101325*400/14.7)
# manifold = Node("Nitrogen", 0.01, 0.1, 287, "manifold")
# fluid_system = Connection(0.0000063, 0, 0)
# network = Network({reg: (copv, manifold), fluid_system: (manifold, amb_node)})
# network.sim(60, 1)
# network.plot_nodes_overlay((copv, manifold), units="E")

# DARCY SPACE FILL VALIDATION
# target mass: 108.86 kg
# n = 24 # num bottles
# vehicle_tank = Node("N2O", 1, 106.5, 293, "Vehicle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, 293, "Fill Tanks")
# amb_node = Ambient()
# fill_line = Connection(0.000006, 0, 0)
# vent_line = Connection(0.0000003, 0, 1, False)
# omv = Connection(0.000008, 0, 0, False)
# test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), omv: (vehicle_tank, amb_node)})
# test_network.sim(1000, 1, {400: (vent_line, True)})
# test_network.plot_nodes_overlay((fill_tanks, vehicle_tank), title="Darcy Space Validation", units="E")

# GN2 TEST
# gas_bottle = Node("Nitrogen", 8, 27, 293)
# amb_node = Ambient()
# # vent_diameter = 3 # mm
# # vent_CdA = 0.8 * vent_diameter**2 * 0.25 * math.pi # mm^2
# vent_CdA = 1.4 * 1e-6
# gas_vent = Connection(vent_CdA, 0, 0, checking=True)
# gas_network = Network({gas_vent : (gas_bottle, amb_node)})
# gas_network.sim(600, 0.1)
# gas_network.plot_nodes_overlay([gas_bottle], units="E")

# ACTIVE PRESS + FILL + FULL THROTTLE
# fill params
# n = 7 # num bottles
# vol = 50 # L
# bottle_temp = 305 # K
# burn_duration = 60 # s
# fill_duration = 1100 # s
# # vehicle tanks
# copv = Node("Nitrogen", 2, 6.61, 293, "COPV") #subscale copv
# ox_tank = Tank("N2O", 0.926*40*1.1, "Nitrogen", 1, 60, 270, 290, "ox_tank")
# fuel_tank = Tank("ISOPROPANOL", 0.309*40*1.1, "Nitrogen", 0.1, 20, 285, 290, "fu_tank")
# amb_node = Ambient()
# chamber = Ambient(P=101325*350*1.25/14.7)
# # valves
# # mdots: OX Mdot = 0.926 kg/s, FU Mdot = 0.309 kg/s
# bb_ox = BangBang(0.000005, 3.4474E+06, 0) # PLV
# bb_fuel = BangBang(0.000005, 3.4474E+06, 0) # PLV
# tv_ox = ThrottleValve(1, target_mdot=0.926*0.6, normal_state=0.926*0.6)
# tv_fuel = ThrottleValve(1, target_mdot=0.309*0.6, normal_state=0.309*0.6)
# # # network and sim
# test_network = Network({bb_ox: (copv, ox_tank), bb_fuel: (copv, fuel_tank), tv_ox: (ox_tank, chamber), tv_fuel: (fuel_tank, chamber)})
# test_network.sim(burn_duration, 1)
# test_network.plot_nodes_overlay((copv, ox_tank, fuel_tank), title=f"bs", units="E")
# test_network.plot_connections_overlay([tv_ox, bb_ox, bb_fuel], units="E")

# --- SIMULATION SETUP ---

# 1. Unit Conversions and Constants
PSI_TO_PA = 6894.75729
L_TO_M3 = 0.001
MM2_TO_M2 = 1e-6

# # Connections
cda_orifice = 1.02 * MM2_TO_M2  # Reasonable restriction for 4500->350 psi regulation

# #PLV VALIDATER
# copv = Node("Nitrogen", 2.27, 26.67, 293, name="COPV")
# amb = Ambient(fluid="Air", P=14*PSI_TO_PA, T=293.15, name="Ambient")
# plv = Connection(cda_orifice)
# nw = Network({plv: (copv, amb)})
# nw.sim(1000, 0.5, verbose_steps=20)
# nw.plot_nodes_overlay([copv], title="Blowdown Simulation: COPV to LOX Tank", units="E")
# nw.plot_connections_overlay([plv], title="Blowdown Simulation: COPV to LOX Tank", units="E")
# Nodes
p_ambient = 300 * PSI_TO_PA
copv = Node("Nitrogen", 4.8, 15.47, 293, name="COPV")
ox_tank = Tank(V_total_L=20, 
            fluid_liq='Oxygen', 
            m_liq=19.43, 
            T_liq=90, 
            fluid_ullage="Nitrogen", 
            P_ullage=500*PSI_TO_PA, 
            T_ullage=150,
            radius=0.2, 
            name="OxTank", htc=500)
fu_tank = Tank(V_total_L=16, 
            fluid_liq='n-Dodecane', 
            m_liq=9.7, 
            T_liq=293, 
            fluid_ullage="Nitrogen", 
            P_ullage=500*PSI_TO_PA, 
            T_ullage=293,
            radius=0.1, 
            name="FuTank", htc=0)
amb = Ambient(fluid="Air", P=300*1.2*PSI_TO_PA, T=293.15, name="Ambient")


# 5. Instantiate Connections
# Orifice: COPV -> Tank
# Connects to the Ullage of the tank (location=1.0)
obb = BangBang(CdA=cda_orifice, 
                    target_pressure=(500*PSI_TO_PA),
                    hysteresis=(5*PSI_TO_PA),
                    location=1.0, 
                    name="Ox Bang-Bang")
fbb = BangBang(CdA=cda_orifice, 
                    target_pressure=(500*PSI_TO_PA),
                    hysteresis=(5*PSI_TO_PA),
                    location=1.0, 
                    name="Fu Bang-Bang")

# Outlet: Tank -> Ambient
# Connects to the Liquid of the tank (location=0.0)
otv = ThrottleValve(1,
                    location=0.0, 
                    name="OxThrottle", normal_state=0.777)
ftv = ThrottleValve(1,
                    location=0.0, 
                    name="FuThrottle", normal_state=0.388)

# throttle profile
actions = {
    0.0: (otv, 0.36),
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
    7.5: (otv, 0.29),
    7.6: (ftv, 0.14),
    8.5: (otv, 0.29),
    8.6: (ftv, 0.14),
    9.5: (otv, 0.29),
    9.6: (ftv, 0.14),
    10.5: (otv, 0.31),
    10.6: (ftv, 0.16),
    11.5: (otv, 0.33),
    11.6: (ftv, 0.16),
    12.5: (otv, 0.33),
    12.6: (ftv, 0.17),
    13.5: (otv, 0.34),
    13.6: (ftv, 0.17),
    14.5: (otv, 0.34),
    14.6: (ftv, 0.17),
    15.5: (otv, 0.29),
    15.6: (ftv, 0.14),
    16.3: (otv, 0.29),
    16.4: (ftv, 0.14),
    17.3: (otv, 0.29),
    17.4: (ftv, 0.14),
    18.3: (otv, 0.29),
    18.4: (ftv, 0.14),
    19.3: (otv, 0.29),
    19.4: (ftv, 0.14),
    20.3: (otv, 0.29),
    20.4: (ftv, 0.14),
    21.3: (otv, 0.29),
    21.4: (ftv, 0.14),
    21.8: (otv, 0.35),
    21.9: (ftv, 0.17),
    22.4: (otv, 0.39),
    22.5: (ftv, 0.2),
    23.1: (otv, 0.39),
    23.2: (ftv, 0.2),
    24.1: (otv, 0.38),
    24.2: (ftv, 0.19),
    25.1: (otv, 0.36),
    25.2: (ftv, 0.18),
    26.1: (otv, 0.34),
    26.2: (ftv, 0.17),
    27.1: (otv, 0.33),
    27.2: (ftv, 0.17),
    28.1: (otv, 0.32),
    28.2: (ftv, 0.16),
    29.1: (otv, 0.32),
    29.2: (ftv, 0.16),
    29.7: (otv, 0.32),
    29.8: (ftv, 0.16)
}
# 6. Define Network Graph
# {Connection: (Upstream_Node, Downstream_Node)}
# Note: Flow direction is automatic based on pressure, but we define topology here.
graph = {
    obb: (copv, ox_tank),
    otv:  (ox_tank, amb),
    fbb: (copv, fu_tank),
    ftv: (fu_tank, amb)
}

network = Network(graph)


# 7. Run Simulation
print("Starting Simulation...")
print(f"Initial State: COPV={copv.P/PSI_TO_PA:.0f} psi, Tank={ox_tank.P/PSI_TO_PA:.0f} psi, Amb={amb.P/PSI_TO_PA:.0f} psi")

dt = 0.1 # 100ms timestep
runtime = 25  # Run for 20 seconds
network.sim(runtime, dt, actions, verbose_steps=10)


# 8. Plotting
# Filter nodes for plotting
plot_nodes = [copv, ox_tank, ox_tank.ullage, fu_tank, fu_tank.ullage]
plot_conns = [obb, otv, fbb, ftv]
network.plot_nodes_overlay(plot_nodes, title="Blowdown Simulation: COPV to LOX Tank", units="E")
network.plot_connections_overlay(plot_conns, title="Blowdown Simulation: COPV to LOX Tank", units="E")

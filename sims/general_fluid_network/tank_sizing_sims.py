from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network, Tank, BangBang, Regulator
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
# burn_duration = 60 # s
# fill_duration = 1100 # s
# vehicle_tank = Node("N2O", 1, vol, 293, "Vehicle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, bottle_temp, "Fill Tanks")
# amb_node = Ambient()
# chamber = Ambient(P=101325*350*1.25/14.7)
# fill_line = Connection(0.00005, 0, 0)
# vent_line = Connection(0.0000005, 0, 1, False)
# tv1 = ThrottleValve(1, target_mdot=0.0, normal_state=0.0)
# test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), tv1: (vehicle_tank, chamber)})
# test_network.sim(fill_duration+burn_duration, 1, {300: (vent_line, True), fill_duration-10: (fill_line, False), fill_duration-1: (vent_line, False), fill_duration: (tv1, 0.57*0.5), fill_duration+10: (tv1, 0.57)})
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
# gas_bottle = Node("Nitrogen", 0.15, 3, 293)
# amb_node = Ambient()
# vent_diameter = 3 # mm
# vent_CdA = 0.8 * vent_diameter**2 * 0.25 * math.pi # mm^2
# vent_CdA = vent_CdA * 1e-6
# gas_vent = Connection(vent_CdA, 0, 0, checking=True)
# gas_network = Network({gas_vent : (gas_bottle, amb_node)})
# gas_network.sim(60, 0.1)
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
# PSI_TO_PA = 6894.75729
# L_TO_M3 = 0.001
# MM2_TO_M2 = 1e-6

# # Connections
# cda_orifice = 5.0 * MM2_TO_M2  # Reasonable restriction for 4500->350 psi regulation

# # Nodes
# p_ambient = 300 * PSI_TO_PA
# copv = Node("Nitrogen", 8.27*9/26.67, 9, 293, name="COPV")
# ox_tank = Tank(V_total_L=20, 
#             fluid_liq='Oxygen', 
#             m_liq=19.43, 
#             T_liq=90, 
#             fluid_ullage="Nitrogen", 
#             P_ullage=500*PSI_TO_PA, 
#             T_ullage=293,
#             radius=0.1, 
#             name="OxTank", htc=200)
# fu_tank = Tank(V_total_L=15.2, 
#             fluid_liq='n-Dodecane', 
#             m_liq=9.7, 
#             T_liq=293, 
#             fluid_ullage="Nitrogen", 
#             P_ullage=500*PSI_TO_PA, 
#             T_ullage=293,
#             radius=0.1, 
#             name="FuTank", htc=0)
# amb = Ambient(fluid="Air", P=300*1.2*PSI_TO_PA, T=293.15, name="Ambient")


# # 5. Instantiate Connections
# # Orifice: COPV -> Tank
# # Connects to the Ullage of the tank (location=1.0)
# obb = BangBang(CdA=cda_orifice, 
#                     target_pressure=(500*PSI_TO_PA),
#                     hysteresis=(5*PSI_TO_PA),
#                     location=1.0, 
#                     name="Ox Bang-Bang")
# fbb = BangBang(CdA=cda_orifice, 
#                     target_pressure=(500*PSI_TO_PA),
#                     hysteresis=(5*PSI_TO_PA),
#                     location=1.0, 
#                     name="Fu Bang-Bang")

# # Outlet: Tank -> Ambient
# # Connects to the Liquid of the tank (location=0.0)
# otv = ThrottleValve(1,
#                     location=0.0, 
#                     name="OxThrottle", normal_state=0.777, target_mdot=0.777)
# ftv = ThrottleValve(1,
#                     location=0.0, 
#                     name="FuThrottle", normal_state=0.388, target_mdot=0.388)


# # 6. Define Network Graph
# # {Connection: (Upstream_Node, Downstream_Node)}
# # Note: Flow direction is automatic based on pressure, but we define topology here.
# graph = {
#     obb: (copv, ox_tank),
#     otv:  (ox_tank, amb),
#     fbb: (copv, fu_tank),
#     ftv: (fu_tank, amb)
# }

# network = Network(graph)


# # 7. Run Simulation
# print("Starting Simulation...")
# print(f"Initial State: COPV={copv.P/PSI_TO_PA:.0f} psi, Tank={ox_tank.P/PSI_TO_PA:.0f} psi, Amb={amb.P/PSI_TO_PA:.0f} psi")

# dt = 0.1 # 100ms timestep
# runtime = 25.0 # Run for 20 seconds
# network.sim(runtime, dt, verbose_steps=250)


# # 8. Plotting
# # Filter nodes for plotting
# plot_nodes = [copv, ox_tank, ox_tank.ullage, fu_tank, fu_tank.ullage]
# network.plot_nodes_overlay(plot_nodes, title="Blowdown Simulation: COPV to LOX Tank", units="E")

# plot_conns = [obb, otv, fbb, ftv]
# network.plot_connections_overlay(plot_conns, title="Connection Flow Rates", units="E")
from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network
#################### TEST CONFIGS #########################
import math


# GAS PULL
# test_node1 = Node("N2O", 35, 60, 293, "Liquid Pull")
# test_node2 = Node("N2O", 35, 60, 293, "Gas Pull")
# amb_node = Ambient()
# test_connection1 = Connection(0.0000063, 0, 0)
# test_connection2 = Connection(0.00001521, 0, 1)
# test_connection3 = Connection(0.0000001742, 0, 0)
# test_network = Network({test_connection1: (test_node1, amb_node), test_connection2: (test_node2, amb_node)})
# test_network.sim(60, 1)
# test_network.plot_nodes_overlay((test_node1, test_node2), units="E")

# FILL + FULL THROTTLE
# target mass: 40 kg
n = 7 # num bottles
vol = 50 # L
bottle_temp = 305 # K
burn_duration = 60 # s
fill_duration = 1100 # s
vehicle_tank = Node("N2O", 1, vol, 293, "Vehicle Tank")
fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, bottle_temp, "Fill Tanks")
amb_node = Ambient()
chamber = Ambient(P=101325*350*1.25/14.7)
fill_line = Connection(0.00005, 0, 0)
vent_line = Connection(0.0000005, 0, 1, False)
tv1 = ThrottleValve(1, target_mdot=0.0, normal_state=0.0)
test_network = Network({fill_line: (fill_tanks, vehicle_tank), vent_line: (vehicle_tank, amb_node), tv1: (vehicle_tank, chamber)})
test_network.sim(fill_duration+burn_duration, 1, {300: (vent_line, True), fill_duration-10: (fill_line, False), fill_duration-1: (vent_line, False), fill_duration: (tv1, 0.57*0.5), fill_duration+10: (tv1, 0.57)})
test_network.plot_nodes_overlay((fill_tanks, vehicle_tank), title=f"{burn_duration}s, {vol}L, {n} Bottles, {bottle_temp}K", units="E")
test_network.plot_connections_overlay([tv1], units="E")

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
# fluid_system = Connection(0.0000063, 0, 0)
# network = Network({reg: (copv, vehicle_tank), fluid_system: (vehicle_tank, amb_node)})
# network.sim(60, 1)
# network.plot_nodes_overlay((copv, vehicle_tank), units="E")

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


# gas_bottle = Node("Nitrogen", 0.15, 3, 293)
# amb_node = Ambient()
# vent_diameter = 3 # mm
# vent_CdA = 0.8 * vent_diameter**2 * 0.25 * math.pi # mm^2
# # vent_CdA = vent_CdA * 1e-6
# gas_vent = Connection(vent_CdA, 0, 0, checking=False)
# gas_network = Network({gas_vent : {gas_bottle, amb_node}})
# gas_network.sim(60, 1)
# gas_network.plot_nodes_overlay([gas_bottle], units="E")


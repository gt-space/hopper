from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network, Regulator, Tank, BangBang, Series, Line, Engine
#################### TEST CONFIGS #########################
import numpy as np
import math
import time
from rocketcea.cea_obj import CEA_Obj
from scipy.optimize import brentq

PSI_TO_PA = 6894.75729
L_TO_M3 = 0.001
MM2_TO_M2 = 1e-6


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
# fill_duration = 1100 # s
# vehicle_tank = Node("N2O", 1, vol, 293, "Vehicle Tank")
# fill_tanks = Node("N2O", 9.0718 * n, 13.4 * n, bottle_temp, "Fill Tanks")
# amb_node = Ambient()
# fill_line = Connection(0.00005, 0, 0)
# vent_line = Connection(0.0000005, 0, 1, False)
# test_network = Network({fill_line: (fill_tanks, amb_node)})
# test_network.sim(fill_duration, 0.01, {150: [(fill_line, True)]})
# test_network.plot_nodes_overlay([fill_tanks])
# test_network.plot_connections_overlay([fill_line], units="E")

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
# test_network.plot_nodes_overlay([vt2], units="E")


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

# # GN2 TEST
# gas_bottle = Node("Nitrogen", 8, 27, 293)
# amb_node = Ambient()
# # vent_diameter = 3 # mm
# # vent_CdA = 0.8 * vent_diameter**2 * 0.25 * math.pi # mm^2
# vent_CdA = 1.4 * 1e-6
# gas_vent = Connection(vent_CdA, 0, 0, checking=True)
# gas_network = Network({gas_vent : (gas_bottle, amb_node)})
# gas_network.sim(600, 0.1)
# gas_network.plot_nodes_overlay([gas_bottle], units="E")

# GN2 FILL
# gas_bottle = Node("Nitrogen", 18, 50, 293)
# flight_tank = Node("Nitrogen", 0.05, 40, 293)
# print(gas_bottle.P/PSI_TO_PA, flight_tank.P/PSI_TO_PA)
# fill_diameter = 0.3 # mm
# fill_CdA = 0.61 * fill_diameter**2 * 0.25 * math.pi # mm^2
# fill_CdA *= 1e-6
# gas_fill = Connection(fill_CdA, 0, 0, checking=True)
# gas_network = Network({gas_fill : (gas_bottle, flight_tank)})
# gas_network.sim(600, 0.1)
# gas_network.plot_nodes_overlay([gas_bottle, flight_tank], units="E")

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

class ThrustProfileTargeter:
    def __init__(self, fuel, oxidizer, MR, At, Ae, Pa, eta_cstar=0.92, eta_cf=0.98):
        self.MR = MR
        self.eta_cstar = eta_cstar
        self.eta_cf = eta_cf
        self.At = At
        self.Ae = Ae
        self.Pa = Pa
        self.eps = Ae / At
        self.cea = CEA_Obj(oxName=oxidizer, fuelName=fuel)

    def get_required_mdots(self, target_thrust_N):
        if target_thrust_N <= 0: return 0.0, 0.0

        def thrust_error(Pc_Pa):
            Pc_psia = max(Pc_Pa / 6894.75729, 14.7)
            
            # Evaluate dynamic efficiencies
            curr_eta_cstar = self.eta_cstar(Pc_psia) if callable(self.eta_cstar) else self.eta_cstar
            curr_eta_cf = self.eta_cf(Pc_psia) if callable(self.eta_cf) else self.eta_cf

            cstar = self.cea.get_Cstar(Pc=Pc_psia, MR=self.MR) * 0.3048 * curr_eta_cstar
            mdot_total = (Pc_Pa * self.At) / cstar

            isp_amb = self.cea.estimate_Ambient_Isp(Pc=Pc_psia, MR=self.MR, eps=self.eps, Pamb=(self.Pa/6894.75729))[0]
            isp = isp_amb * curr_eta_cstar * curr_eta_cf
            
            return (mdot_total * isp * 9.81) - target_thrust_N

        try:
            Pc_solution = brentq(thrust_error, self.Pa * 1.01, 200 * 101325.0)
            Pc_psia = Pc_solution / 6894.75729
            
            curr_eta_cstar = self.eta_cstar(Pc_psia) if callable(self.eta_cstar) else self.eta_cstar
            cstar = self.cea.get_Cstar(Pc=Pc_psia, MR=self.MR) * 0.3048 * curr_eta_cstar
            mdot_total = (Pc_solution * self.At) / cstar
            
            mdot_fu = mdot_total / (1.0 + self.MR)
            return mdot_total - mdot_fu, mdot_fu
        except ValueError:
            return 0.0, 0.0

    def build_action_dict(self, thrust_profile, otv_obj, ftv_obj, otv_max_mdot, ftv_max_mdot):
        actions = {}
        for t, thrust in thrust_profile.items():
            m_ox, m_fu = self.get_required_mdots(thrust)
            actions[t] = [
                (otv_obj, round(m_ox / otv_max_mdot, 3)),
                (ftv_obj, round(m_fu / ftv_max_mdot, 3))
            ]
        return actions
        
    
# # # Connections
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
copv = Node("Nitrogen", 5.1, 18, 293, name="COPV")
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
obb = BangBang(CdA=2*cda_orifice, 
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
# otv = ThrottleValve(10e-6,
#                     location=0.0, 
#                     name="OxThrottle", normal_state=1)
# ftv = ThrottleValve(6e-6,
#                     location=0.0, 
#                     name="FuThrottle", normal_state=1)

# The first argument is now max_mdot (kg/s) because mode="mdot"
otv = ThrottleValve(0.777, 
                    location=0.0, 
                    normal_state=1,
                    mode="mdot",
                    name="OxThrottle")

ftv = ThrottleValve(0.388, 
                    location=0.0, 
                    normal_state=1,
                    mode="mdot",
                    name="FuThrottle")
# injectors:
# ox_inj = Connection(32.3e-6)
# fu_inj = Connection(11.2e-6)

# Lines:
obb_line = Line(0.00492, 1, 1.5e-5, location=1)
fbb_line = Line(0.00492, 1, 1.5e-5, location=1)
otv_line = Line(0.009, 1, 2e-4, location=0)
ftv_line = Line(0.009, 0.1, 2e-4, location=0)

obb_series = Series([obb_line, obb], "obb_series")
fbb_series = Series([fbb_line, fbb], "fbb_series")
otv_series = Series([otv_line, otv], "otv_series")
ftv_series = Series([ftv_line, ftv], "ftv_series")


# --- Define Efficiency Curves ---
def dynamic_eta_cstar(Pc_psi):
    # Example: 92% above 250 psi, dropping to 75% at 100 psi
    if Pc_psi >= 250: return 0.92
    return 0.75 + (0.92 - 0.75) * max(0, (Pc_psi - 100) / 150)

def dynamic_eta_cf(Pc_psi):
    # Example: 96% above 250 psi, dropping sharply due to over-expansion
    if Pc_psi >= 250: return 0.96
    return 0.80 + (0.96 - 0.80) * max(0, (Pc_psi - 100) / 150)

# 1. Define your engine properties
targeter = ThrustProfileTargeter(
    fuel="Kerosene", 
    oxidizer="LOX", 
    MR=1.8,                 
    eta_cstar=dynamic_eta_cstar, # Passing the function 
    eta_cf=dynamic_eta_cf,       # Passing the function
    At=0.00041968,  
    Ae=0.00268451,  
    Pa=14.8 * PSI_TO_PA, 
)

# 2. Design your Thrust Curve (Time : Thrust in Newtons)
thrust_curve = {
    0.0: 1115.9,
    0.1: 1176.9,
    0.3: 1238.0,
    0.6: 1309.1,
    1.0: 1365.7,
    1.6: 1382.5,
    2.5: 1347.9,
    3.5: 1276.7,
    4.5: 1209.9,
    5.5: 1007.4,
    6.5: 890.0,
    7.5: 890.0,
    8.5: 890.0,
    9.5: 890.0,
    10.5: 970.4,
    11.5: 1015.4,
    12.5: 1034.7,
    13.5: 1040.9,
    14.5: 1040.6,
    15.5: 896.9,
    16.3: 890.0,
    17.3: 890.0,
    18.3: 890.0,
    19.3: 890.0,
    20.3: 890.0,
    21.3: 890.0,
    21.8: 1074.1,
    22.4: 1213.5,
    23.1: 1218.3,
    24.1: 1172.4,
    25.1: 1115.0,
    26.1: 1066.2,
    27.1: 1030.0,
    28.1: 1004.8,
    29.1: 987.6,
    29.7: 979.9
}
# 3. Automatically generate the exact mdot actions needed
actions = targeter.build_action_dict(thrust_curve, otv, ftv, 0.777, 0.388)

# ... (network graph definition remains the same) ...

engine = Engine(fuel="Kerosene", 
                oxidizer="LOX", 
                ox_conn=otv_series,  
                fuel_conn=ftv_series, 
                eta_cstar=dynamic_eta_cstar, # Passing the function
                eta_cf=dynamic_eta_cf,       # Passing the function
                At=0.00041968,  
                Ae=0.00268451,   
                Pa=14.8 * PSI_TO_PA, 
                name="Rocket Engine")    

# throttle profile
# actions = {
#     0.0: [(otv, 0.8)],
#     0.1: [(otv, 0.844), (ftv, 0.864)],
#     0.3: [(otv, 0.889), (ftv, 0.909)],
#     0.6: [(otv, 0.933), (ftv, 0.955)],
#     1.0: [(otv, 0.978), (ftv, 1.0)],
#     1.6: [(otv, 1.0), (ftv, 1.0)],
#     2.5: [(otv, 0.978), (ftv, 1.0)],
#     3.5: [(otv, 0.911), (ftv, 0.955)],
#     4.5: [(otv, 0.867), (ftv, 0.909)],
#     5.5: [(otv, 0.733), (ftv, 0.727)],
#     6.5: [(otv, 1), (ftv, 1)],
#     7.5: [(otv, 1), (ftv, 1)],
#     8.5: [(otv, 1), (ftv, 1)],
#     9.5: [(otv, 0.644), (ftv, 0.636)],
#     10.5: [(otv, 0.689), (ftv, 0.727)],
#     11.5: [(otv, 0.733), (ftv, 0.727)],
#     12.5: [(otv, 0.733), (ftv, 0.773)],
#     13.5: [(otv, 0.756), (ftv, 0.773)],
#     14.5: [(otv, 0.756), (ftv, 0.773)],
#     15.5: [(otv, 0.644), (ftv, 0.636)],
#     16.3: [(otv, 0.644), (ftv, 0.636)],
#     17.3: [(otv, 0.644), (ftv, 0.636)],
#     18.3: [(otv, 0.644), (ftv, 0.636)],
#     19.3: [(otv, 0.644), (ftv, 0.636)],
#     20.3: [(otv, 0.644), (ftv, 0.636)],
#     21.3: [(otv, 0.644), (ftv, 0.636)],
#     21.8: [(otv, 0.778), (ftv, 0.773)],
#     22.4: [(otv, 0.867), (ftv, 0.909)],
#     23.1: [(otv, 0.867), (ftv, 0.909)],
#     24.1: [(otv, 0.844), (ftv, 0.864)],
#     25.1: [(otv, 0.8), (ftv, 0.818)],
#     26.1: [(otv, 0.756), (ftv, 0.773)],
#     27.1: [(otv, 0.733), (ftv, 0.773)],
#     28.1: [(otv, 0.711), (ftv, 0.727)],
#     29.1: [(otv, 0.711), (ftv, 0.727)],
#     29.7: [(otv, 0.711), (ftv, 0.727)]
# }
# 6. Define Network Graph
# {Connection: (Upstream_Node, Downstream_Node)}
# Note: Flow direction is automatic based on pressure, but we define topology here.

graph = {
    obb_series: (copv, ox_tank),
    otv_series:  (ox_tank, engine),
    fbb_series: (copv, fu_tank),
    ftv_series: (fu_tank, engine)
}

network = Network(graph)

# 7. Run Simulation
print("Starting Simulation...")
print(f"Initial State: COPV={copv.P/PSI_TO_PA:.0f} psi, Tank={ox_tank.P/PSI_TO_PA:.0f} psi, Amb={amb.P/PSI_TO_PA:.0f} psi")

dt = 0.1 # 100ms timestep
runtime = 30  # Run for 20 seconds
t = time.time()
network.sim(runtime, dt, actions, verbose_steps=0)
print(t-time.time())

# 8. Plotting
# Filter nodes for plotting
plot_nodes = [copv, ox_tank, ox_tank.ullage, fu_tank, fu_tank.ullage, engine]
plot_conns_ox = [otv_series]
plot_conns_fu = [ftv_series]

network.plot_nodes_overlay(plot_nodes, title="Blowdown Simulation: COPV to LOX Tank", units="E")
network.plot_connections_overlay(plot_conns_ox, title="Blowdown Simulation: COPV to LOX Tank", units="E")
network.plot_connections_overlay(plot_conns_fu, title="Blowdown Simulation: COPV to LOX Tank", units="E")
network.plot_connections_overlay([obb_series, fbb_series], title="Blowdown Simulation: COPV to LOX Tank", units="E")

import matplotlib.pyplot as plt

# Custom plot for Engine Performance
fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
time = engine.history['time']

axs[0].plot(time, np.array(engine.history['P']) / PSI_TO_PA, color='red')
axs[0].set_ylabel("Chamber Pressure [psi]")
axs[0].grid(True)

axs[1].plot(time, engine.history['thrust'], color='orange')
axs[1].set_ylabel("Thrust [N]")
axs[1].grid(True)

axs[2].plot(time, engine.history['MR'], color='purple')
axs[2].set_ylabel("Mixture Ratio")
axs[2].set_xlabel("Time [s]")
axs[2].grid(True)

plt.suptitle("Engine Transient Performance")
plt.show()

# --- Plot Operating Box ---
# Filter out startup/shutdown transients (Thrust < 100 N) to keep the box clean
valid_idx = [i for i, thrust in enumerate(engine.history['thrust']) if thrust > 100]
op_Pc = [engine.history['P'][i] / PSI_TO_PA for i in valid_idx]
op_MR = [engine.history['MR'][i] for i in valid_idx]
op_time = [engine.history['time'][i] for i in valid_idx]

fig, ax = plt.subplots(figsize=(8, 6))

# Define your redline limits (Adjust these to your engine's actual limits)
min_Pc, max_Pc = 150, 400
min_MR, max_MR = 1.2, 2.5

# Draw the Box
ax.axvspan(min_MR, max_MR, ymin=0, ymax=1, color='green', alpha=0.1, label='Safe Operating Box')
ax.axhline(min_Pc, color='red', linestyle='--', linewidth=1)
ax.axhline(max_Pc, color='red', linestyle='--', linewidth=1)
ax.axvline(min_MR, color='red', linestyle='--', linewidth=1)
ax.axvline(max_MR, color='red', linestyle='--', linewidth=1)

# Scatter plot the engine's trajectory colored by time
sc = ax.scatter(op_MR, op_Pc, c=op_time, cmap='viridis', s=10, zorder=5)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Time [s]')

ax.set_xlim(min_MR - 0.2, max_MR + 0.2)
ax.set_ylim(min_Pc - 50, max_Pc + 50)
ax.set_xlabel("Mixture Ratio (O/F)")
ax.set_ylabel("Chamber Pressure [psi]")
ax.set_title("Engine Operating Envelope")
ax.legend(loc='upper left')
ax.grid(True)

plt.show()
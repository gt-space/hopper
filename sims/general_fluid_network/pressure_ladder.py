from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network, Regulator, Tank, BangBang, PropsSI_auto
#################### TEST CONFIGS #########################
import math
print(PropsSI_auto('D', 'T', 300, 'P', 101325, 'Nitrogen'))
# --- SIMULATION SETUP ---
 
# Ref Lines
# Line	Segment	Length (m)
# LOX	Tank → OTV	0.7735 m
# LOX	OTV → TCA	0.1524 m
# Fuel	Tank → FTV	0.3048 m
# Fuel	FTV → TCA	0.3302 m
 
# Fuel Side: AL6061
# Tube OD	0.375	in
# Tube Thickness	0.028	in
 
# Ox Side: Stainless 304
# Tube OD	0.375	in
# Tube Thickness	0.01	in

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
                    name="OxThrottle", normal_state=0.36)
ftv = ThrottleValve(1,
                    location=0.0, 
                    name="FuThrottle", normal_state=0.19)

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

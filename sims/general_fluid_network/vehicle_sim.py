import numpy as np
import time
import csv
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from scipy.optimize import brentq
from general_fluid_network import Node, Ambient, Connection, ThrottleValve, Network, Regulator, Tank, BangBang, Series, Line, Engine

# --- CONSTANTS ---
PSI_TO_PA = 6894.75729
L_TO_M3 = 0.001
MM2_TO_M2 = 1e-6
IN_TO_M = 0.0254

# =============================================================================
# ENGINE EFFICIENCY CURVES AS A FUNCTION OF PC
# =============================================================================
def dynamic_eta_cstar(Pc_psi):
    """Cstar efficiency based on chamber pressure."""
    return 0.8 + (0.85 - 0.8) * max(1, (Pc_psi - 100) / 200)

def dynamic_eta_cf(Pc_psi):
    """Cf efficiency based on chamber pressure."""
    return 0.94 + (0.97 - 0.94) * max(1, (Pc_psi - 100) / 150)

# =============================================================================
# TARGETER CLASS
# =============================================================================
class ThrustProfileTargeter:
    """
    Reverse-calculates exact mdot_ox, mdot_fu, and Pc required to hit a specific 
    thrust profile over time at a fixed Mixture Ratio.
    """
    def __init__(self, fuel, oxidizer, MR, At, Ae, Pa, eta_cstar, eta_cf):
        self.MR = MR
        self.eta_cstar = eta_cstar
        self.eta_cf = eta_cf
        self.At = At
        self.Ae = Ae
        self.Pa = Pa
        self.eps = Ae / At
        self.cea = CEA_Obj(oxName=oxidizer, fuelName=fuel)

    def get_required_mdots(self, target_thrust_N):
        if target_thrust_N <= 0: return 0.0, 0.0, 14.7

        def thrust_error(Pc_Pa):
            Pc_psia = max(Pc_Pa / PSI_TO_PA, 14.7)
            
            curr_eta_cstar = self.eta_cstar(Pc_psia) if callable(self.eta_cstar) else self.eta_cstar
            curr_eta_cf = self.eta_cf(Pc_psia) if callable(self.eta_cf) else self.eta_cf

            cstar = self.cea.get_Cstar(Pc=Pc_psia, MR=self.MR) * 0.3048 * curr_eta_cstar
            mdot_total = (Pc_Pa * self.At) / cstar

            isp_amb = self.cea.estimate_Ambient_Isp(Pc=Pc_psia, MR=self.MR, eps=self.eps, Pamb=(self.Pa/PSI_TO_PA))[0]
            isp = isp_amb * curr_eta_cstar * curr_eta_cf
            
            return (mdot_total * isp * 9.81) - target_thrust_N

        try:
            Pc_solution = brentq(thrust_error, self.Pa * 1.01, 200 * 101325.0)
            Pc_psia = Pc_solution / PSI_TO_PA
            
            curr_eta_cstar = self.eta_cstar(Pc_psia) if callable(self.eta_cstar) else self.eta_cstar
            cstar = self.cea.get_Cstar(Pc=Pc_psia, MR=self.MR) * 0.3048 * curr_eta_cstar
            mdot_total = (Pc_solution * self.At) / cstar
            
            mdot_fu = mdot_total / (1.0 + self.MR)
            mdot_ox = mdot_total - mdot_fu
            return mdot_ox, mdot_fu, Pc_psia
        except ValueError:
            print(f"Solver failed for target thrust: {target_thrust_N} N")
            return 0.0, 0.0, 14.7

    def build_action_dict(self, thrust_profile, otv_obj, ftv_obj, otv_max_mdot, ftv_max_mdot):
        actions = {}
        for t, thrust in thrust_profile.items():
            m_ox, m_fu, _ = self.get_required_mdots(thrust)
            # Removed the round(..., 3) to prevent quantization drifting the MR
            actions[t] = [
                (otv_obj, m_ox / otv_max_mdot),
                (ftv_obj, m_fu / ftv_max_mdot)
            ]
        return actions

    def export_csv(self, thrust_profile, filename="engine_targets.csv"):
        print(f"Exporting thrust targets to {filename}...")
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Time [s]", "Target Thrust [N]", "Target Pc [psi]", "Mdot Ox [kg/s]", "Mdot Fu [kg/s]", "Mdot Total [kg/s]"])
            for t, thrust in thrust_profile.items():
                m_ox, m_fu, pc_psi = self.get_required_mdots(thrust)
                writer.writerow([t, thrust, round(pc_psi, 2), round(m_ox, 4), round(m_fu, 4), round(m_ox+m_fu, 4)])

# =============================================================================
# SIMULATION SETUP
# =============================================================================
# --- 1. Nodes ---
p_ambient = 300 * PSI_TO_PA
amb = Ambient(fluid="Air", P=300*1.2*PSI_TO_PA, T=293.15, name="Ambient")
copv = Node("Nitrogen", 4.8, 15.47, 293, name="COPV")

ox_tank = Tank(V_total_L=20, fluid_liq='Oxygen', m_liq=19.43, T_liq=90, 
                fluid_ullage="Nitrogen", P_ullage=500*PSI_TO_PA, T_ullage=150,
                radius=0.2, name="OxTank", htc=500)
                
fu_tank = Tank(V_total_L=16, fluid_liq='n-Dodecane', m_liq=9.7, T_liq=293, 
                fluid_ullage="Nitrogen", P_ullage=500*PSI_TO_PA, T_ullage=293,
                radius=0.1, name="FuTank", htc=0)

# --- 2. Targeter & Profile ---
targeter = ThrustProfileTargeter(
    fuel="Kerosene", oxidizer="LOX", MR=1.8, 
    eta_cstar=dynamic_eta_cstar, eta_cf=dynamic_eta_cf, 
    At=0.00067226, Ae=0.00235559904, Pa=14.8 * PSI_TO_PA
)

# from 1dof thrust profile
thrust_curve_1dof = {
    0.0: 1115.9, 0.1: 1176.9, 0.3: 1238.0, 0.6: 1309.1, 1.0: 1365.7, 
    1.6: 1382.5, 2.5: 1347.9, 3.5: 1276.7, 4.5: 1209.9, 5.5: 1007.4, 
    6.5: 890.0,  7.5: 890.0,  8.5: 890.0,  9.5: 890.0,  10.5: 970.4, 
    11.5: 1015.4, 12.5: 1034.7, 13.5: 1040.9, 14.5: 1040.6, 15.5: 896.9, 
    16.3: 890.0, 17.3: 890.0, 18.3: 890.0, 19.3: 890.0, 20.3: 890.0, 
    21.3: 890.0, 21.8: 1074.1, 22.4: 1213.5, 23.1: 1218.3, 24.1: 1172.4, 
    25.1: 1115.0, 26.1: 1066.2, 27.1: 1030.0, 28.1: 1004.8, 29.1: 987.6, 29.7: 979.9
}

thrust_curve_throttle_test = {
    0.0: 889.64,
    0.1: 898.53,
    0.3: 916.33,
    0.6: 943.02,
    1.0: 978.61,
    1.6: 1031.99,
    2.5: 1112.06,
    3.5: 1201.02,
    4.5: 1289.98,
    5.5: 1378.95,
    6.5: 1467.91,
    7.5: 1556.88,
    8.5: 1645.84,
    9.5: 1734.81,
    10.5: 1823.77,
    11.5: 1912.73,
    12.5: 2001.70,
    13.5: 2090.66,
    14.5: 2179.63,
    15.5: 2224.11
}

thrust_curve = thrust_curve_throttle_test
end_time = next(reversed(thrust_curve))

# --- 3. Connections ---
# bang bang orifice
cda_orifice = 1.02 * MM2_TO_M2

# Pressurization
obb_line = Line(0.00492, 1, 1.5E-5, location=1)
obb = BangBang(CdA=2*cda_orifice, target_pressure=(500*PSI_TO_PA), hysteresis=(5*PSI_TO_PA), location=1.0, name="Ox Bang-Bang")
obb_series = Series([obb_line, obb], "obb_series")

fbb_line = Line(0.00492, 1, 1.5E-5, location=1)
fbb = BangBang(CdA=cda_orifice, target_pressure=(500*PSI_TO_PA), hysteresis=(5*PSI_TO_PA), location=1.0, name="Fu Bang-Bang")
fbb_series = Series([fbb_line, fbb], "fbb_series")

# Engine Feed (Lines + Valves + Injectors in Series)
otv_max, ftv_max = 0.777, 0.388

# Tubes
ID_3_8_028 = ((3/8) - 2*0.028) * IN_TO_M # ID of 3/8" Tube, WT = 0.028", [in] --> [m]
ID_1_4_028 = ((1/4) - 2*0.028) * IN_TO_M # ID of 1/4" Tube, WT = 0.028" [in] --> [m]


# 1. Tank --> Throttle Valves
otv_line = Line(ID_3_8_028, 1, 1.5E-5, location=0, name="OX Tank to OTV")
ftv_line = Line(ID_1_4_028, 0.3048, 1.5E-5, location=0, name="FU Tank to FTV")

# 2. Throttle Valves
otv = ThrottleValve(otv_max, location=0.0, normal_state=1.0, mode="mdot", name="OTV")
ftv = ThrottleValve(ftv_max, location=0.0, normal_state=1.0, mode="mdot", name="FTV")

# 3. Throttle Valves --> Injectors
otv_to_inj =  Line(ID_3_8_028, 0.2032, 1.5E-5, location=0, name="OX Tank to OINJ")
ftv_to_inj = Line(ID_1_4_028, 0.2032, 1.5E-5, location=0, name="FU Tank to FINJ")

# 4. Injectors 
ox_inj = Connection(15.3E-6, name="OINJ")
fu_inj = Connection(10.6E-6, name="FINJ")

# 5. Regen
fu_regen = Connection(25e-6, name="FREGEN")

# 6. Combine
otv_series = Series([otv_line, otv, otv_to_inj, ox_inj], "otv_series")
ftv_series = Series([ftv_line, ftv, fu_inj, ftv_to_inj, fu_regen], "ftv_series")

# Generate Action Dict
actions = targeter.build_action_dict(thrust_curve, otv, ftv, otv_max, ftv_max)

# --- 4. Engine & Network ---
engine = Engine(fuel="Kerosene", oxidizer="LOX", 
                ox_conn=otv_series, fuel_conn=ftv_series, # Safely use the Series blocks now
                eta_cstar=dynamic_eta_cstar, eta_cf=dynamic_eta_cf, 
                At=0.00067226, Ae=0.00235559904, Pa=14.8 * PSI_TO_PA, 
                name="Rocket Engine")

graph = {
    obb_series: (copv, ox_tank),
    fbb_series: (copv, fu_tank),
    otv_series: (ox_tank, engine),
    ftv_series: (fu_tank, engine)
}

network = Network(graph)

# --- 5. Run Sim ---
print("Starting Simulation...")
dt = 0.1 
runtime = end_time  
t0 = time.time()
network.sim(runtime, dt, actions, verbose_steps=0)
print(f"Simulation complete in {time.time()-t0:.2f} seconds.")

# =============================================================================
# EXPORT ACTUAL PERFORMANCE TO CSV
# =============================================================================
export_filename = "actual_performance.csv"
print(f"Exporting actual simulated performance to {export_filename}...")
with open(export_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write CSV Headers
    writer.writerow([
        "Time [s]", 
        "Actual Pc [psi]", 
        "Actual Thrust [N]", 
        "Actual MR", 
        "Mdot Ox [kg/s]", 
        "Mdot Fu [kg/s]", 
        "Total Mdot [kg/s]",
        "Actual C* [m/s]", 
        "Actual Isp [s]"
    ])
    
    # Loop through the engine history and write each row
    for i in range(len(engine.history['time'])):
        t_val = engine.history['time'][i]
        pc_val = engine.history['P'][i] / PSI_TO_PA
        thrust_val = engine.history['thrust'][i]
        mr_val = engine.history['MR'][i]
        m_ox_val = engine.history['mdot_ox'][i]
        m_fu_val = engine.history['mdot_fu'][i]
        cstar_val = engine.history['cstar'][i]
        isp_val = engine.history['Isp'][i]
        
        writer.writerow([
            round(t_val, 3),
            round(pc_val, 2),
            round(thrust_val, 2),
            round(mr_val, 4),
            round(m_ox_val, 4),
            round(m_fu_val, 4),
            round(m_ox_val + m_fu_val, 4),
            round(cstar_val, 2),
            round(isp_val, 2)
        ])
        
print("Export complete.")

# Post-Processing
ox_stiffness = (np.array(ox_inj.history['dP']) / np.array(engine.history['P'])) * 100
fu_stiffness = (np.array(fu_inj.history['dP']) / np.array(engine.history['P'])) * 100

# =============================================================================
# PLOTTING
# =============================================================================

# --- 1. Fluid Network Plots ---
plot_nodes = [copv, ox_tank, ox_tank.ullage, fu_tank, fu_tank.ullage, engine]
plot_conns_ox = [otv_series]
plot_conns_fu = [ftv_series]

network.plot_nodes_overlay(plot_nodes, title="Network Simulation: Nodes", units="E")
network.plot_connections_overlay(plot_conns_ox, title="Network Simulation: Ox Feed System", units="E")
network.plot_connections_overlay(plot_conns_fu, title="Network Simulation: Fuel Feed System", units="E")
network.plot_connections_overlay([obb_series, fbb_series], title="Network Simulation: Pressurization System", units="E")

# --- 2. Transient Engine Performance ---
fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
time_arr = engine.history['time']

axs[0].plot(time_arr, np.array(engine.history['P']) / PSI_TO_PA, color='red')
axs[0].set_ylabel("Chamber Pressure [psi]")
axs[0].grid(True)

axs[1].plot(time_arr, engine.history['thrust'], color='orange', label="Simulated Thrust")
# Overlay the target curve for verification
axs[1].plot(list(thrust_curve.keys()), list(thrust_curve.values()), color='black', linestyle='--', label="Target Thrust")
axs[1].set_ylabel("Thrust [N]")
axs[1].legend()
axs[1].grid(True)
axs[2].plot(time_arr, engine.history['MR'], color='purple')
axs[2].set_ylabel("Mixture Ratio")
axs[2].set_xlabel("Time [s]")
axs[2].set_ylim(1.0, 2.5) 
axs[2].grid(True)

plt.suptitle("Engine Transient Performance")
plt.tight_layout()

# Stifness Plots
plt.figure()
plt.plot(time_arr, ox_stiffness, color='g', label='OINJ Stiffness')
plt.plot(time_arr, fu_stiffness, color='r', label='FINJ Stiffness')
plt.legend()
plt.xlabel("Time (sec)")
plt.ylabel("Stiffness (%)")
plt.title("Injector Stiffness")

plt.show() # Display engine performance plot
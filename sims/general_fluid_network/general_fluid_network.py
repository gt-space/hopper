import os
import numpy as np
import matplotlib.pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import CoolProp.CoolProp as CP
from rocketcea.cea_obj import CEA_Obj
from scipy.optimize import fsolve

# ==============================================================================
# EOS SOLVER INITIALIZATION
# ==============================================================================

try:
    RP = REFPROPFunctionLibrary('C:\\Program Files\\REFPROP') # Modify your install location if necessary
    RP.SETPATHdll(os.environ.get('RPPREFIX', r"C:\Program Files\REFPROP")) # Modify your install location if necessary
    REFPROP = True
    print("Using REFPROP.")
except ValueError:
    pass
    REFPROP = False
    print("Using CoolProp.")
    REFPROP=False

def PropsSI_auto(output: str, key1: str, val1: float, key2: str, val2: float, fluid: str):
    """
    Selects a fluid EOS solver depending if you have a REFPROP license. Otherwise, CoolProp will be used.
    """
    if REFPROP:
        if output == "Q":
            result = RP.REFPROPdll(
            fluid,
            key1 + key2,
            "QMASS",
            RP.MASS_BASE_SI,  # base SI units
            0, 0,             # iFlag, iUnits
            val1, val2,
            [1.0]              # composition (pure fluid)
            )
        else:
            result = RP.REFPROPdll(
            fluid,
            key1 + key2,
            output,
            RP.MASS_BASE_SI,  # base SI units
            0, 0,             # iFlag, iUnits
            val1, val2,
            [1.0]              # composition (pure fluid)
            )

        return result.Output[0]
    else:
        if output == "E":
            output = "U"
        elif output == "CV":
            output = "CVMASS"
        elif output == "CP":
            output = "CPMASS"
        if key1 == "E":
            key1 = "U"
        if key2 == "E":
            key2 = "U"
        return  CP.PropsSI(output, key1, val1, key2, val2, fluid)

# ==============================================================================
# AUXILLARY CLASSES
# ==============================================================================
"""
TODO:
CREATE HEAT TRANSFER CLASS: IDEA IS THAT IT CAN BE ADDED TO A NODE OR CONNECTION AND CAN EFFECTIVELY CREATE A THERMAL NETWORK
ACTUALLY SHOULD JUST ADD A THERMAL NETWORK THAT YOU CAN JUST STICK TO NODES

CREATE CONTROLLER CLASS AND SUBCLASSES:
CAN ATTACH TO VALVE-LIKE CONNECTIONS 
UPDATE IN NETWORK SIM TO CONTROL VALVE OPENING/THROTTLE BASED ON WHATEVER CONTROL SCHEME
BANG BANG
THROTTLE
PID CONTROL
REDLINES
"""
class HeatTranfer():
    pass


class Controller():
    def __init__(self):
        pass
    

# ==============================================================================
# NODE CLASS AND SUBCLASSES
# ==============================================================================

class Node():
    """
    Node class. State defined by total density d (kg/m^3) and enthalpy K(J).
    Initialized by fluid, mass m (kg), volume V (L), tempurature T (K), and name.
    """
    def __init__(self, fluid, m, V, T, name="node", type="m"):
        if type == "m":
            self.fluid = fluid
            self.m = float(m) # node mass [kg]
            self.V = float(V) / 1000.0  # convert L -> m^3
            self.name = name

            # Initialize state using T and density computed from m/V
            self.d = self.m / self.V
            # specific enthalpy from (D,T)
            u_spec = PropsSI_auto('E', 'D', self.d, 'T', float(T), self.fluid)
            self.U = u_spec * self.m
            # derived (will also populate m_l, m_v)
            self._flash_from_DU(self.d, self.U)
            self.H = self.h * self.m

            # node history dict initialization
            self.history = {k: [] for k in ["time","Q","P","T","U","h","d","m","m_l","m_v", "fill_level", "s"]}

    def _flash_from_DU(self, d, U):
        """
        Given bulk density d (kg/m3) and total enthalpy H (J),
        compute T, P, h and split m into m_l, m_v if two-phase.
        """
        m = self.m
        if m <= 0:
            # safety floor
            m = 1e-12
            self.m = m

        u = U / m  # specific enthalpy J/kg
        # try to get T and P from (D,H)
        try:
            T = PropsSI_auto('T', 'D', d, 'E', u, self.fluid)
            P = PropsSI_auto('P', 'D', d, 'E', u, self.fluid)
            s = PropsSI_auto('S', 'D', d, 'E', u, self.fluid)
            h = PropsSI_auto('H', 'D', d, 'E', u, self.fluid)
            self.phase = CP.PhaseSI('D', d, 'H', h, self.fluid) # use only CoolProp here, REFPROP phase lookup behaves weirdly
        except Exception as e:
            # Fallback for extreme states or errors
            # Attempt to recover using P-H or T-Q if possible, otherwise raise
            raise RuntimeError(f"CoolProp lookup failed in flash: d={d}, h={u}, err={e}") from e

        self.T = T
        self.P = P
        self.u = u
        self.d = d
        self.s = s
        self.h = h

        if self.phase == "twophase":
            Q = PropsSI_auto('Q', 'D', d, 'E', u, self.fluid)  # 0-1
            # saturated liquid and vapor properties at P
            h_l = PropsSI_auto('H', 'P', P, 'Q', 0, self.fluid)
            h_v = PropsSI_auto('H', 'P', P, 'Q', 1, self.fluid)
            d_l = PropsSI_auto('D', 'P', P, 'Q', 0, self.fluid)
            d_v = PropsSI_auto('D', 'P', P, 'Q', 1, self.fluid)

            self.Q = Q
            self.h_l = h_l
            self.h_v = h_v
            self.d_l = d_l
            self.d_v = d_v

            # masses
            self.m_v = Q * self.m
            self.m_l = self.m - self.m_v

            # fill level (volume fraction of liquid in tank)
            # liquid volume = m_l / d_l
            self.fill_level = (self.m_l / self.d_l) / self.V
        else:
            # single phase (liquid or gas)
            self.Q = None
            self.m_v = self.m if self.phase in ("gas", "supercritical", "supercritical_gas") else 0.0
            self.m_l = self.m - self.m_v

            # set phase-specific properties equal to bulk
            self.h_l = self.h_v = self.h
            self.d_l = self.d_v = self.d
            self.fill_level = 1.0 if self.phase in ("liquid", "supercritical_liquid") else 0.0

        # safe Cp/Cv/gamma/R in single-phase gas
        try:
            self.Cp = PropsSI_auto('CP', 'D', self.d, 'H', self.h, self.fluid)
            self.Cv = PropsSI_auto('CV', 'D', self.d, 'H', self.h, self.fluid)
            self.gamma = self.Cp / self.Cv if (self.Cv and self.Cp) else None
            self.R = self.Cp - self.Cv if (self.Cp and self.Cv) else None
        except Exception:
            self.Cp = self.Cv = self.gamma = self.R = None

    def update(self, mdot, Hdot, dt):
        """
        Updates node state based on an input mdot (kg/s), an input Hdot (J/s),
        as well as the sim timestep dt (s).
        """
        # apply updates
        self.m += mdot * dt
        self.U += Hdot * dt

        # numerical safety
        if self.m < 1e-12:
            self.m = 1e-12

        # recompute density and flash to get phase split
        d_new = self.m / self.V
        self._flash_from_DU(d_new, self.U)

    def log_state(self, t=0.0):
        """
        Log node state at each timestep throughout a network sim.
        """
        self.history["time"].append(t)
        self.history["Q"].append(self.Q)
        self.history["P"].append(self.P)
        self.history["T"].append(self.T)
        self.history["U"].append(self.U)
        self.history["h"].append(self.h)
        self.history["d"].append(self.d)
        self.history["m"].append(self.m)
        self.history["m_l"].append(self.m_l)
        self.history["m_v"].append(self.m_v)
        self.history["fill_level"].append(self.fill_level)
        self.history["s"].append(self.s)


class Ambient(Node):
    """
    Subclass of Node to represnt ambient properties. 
    Unchanging regardless of updates into or out of it.
    """
    def __init__(self, fluid="Air", P=101325, T=293.15, name="ambient"):
        d = PropsSI_auto("D", "P", P, "T", T, fluid)
        super().__init__(fluid, m=1.0, V=1000.0/d, T=T, name=name)

    def update(self, mdot, Hdot, dt):
        """
        Ignore mass/energy inflows, hold fixed at initial state.
        """        
        pass


class Manifold(Node):
    """
    Subclass of Node to represent a volumeless manifold.
    """
    # TODO
    pass


class Engine:
    """
    Dynamic Engine Node-ish model.
    Calculates Pc quasi-steadily based on incoming mdot from specified ox and fuel connections.
    Outputs performance parameters like Thrust and Isp dynamically using RocketCEA.
    """
    def __init__(self, fuel, oxidizer, ox_conn, fuel_conn, At, Ae, Pa, eta_cstar=0.92, eta_cf=0.98, name="engine"):
        self.name = name
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.ox_conn = ox_conn
        self.fuel_conn = fuel_conn
        self.eta_cstar = eta_cstar # Can be float or callable: f(Pc_psia)
        self.eta_cf = eta_cf       # Can be float or callable: f(Pc_psia)
        self.At = At
        self.Ae = Ae
        self.Pa = Pa
        
        self.cea = CEA_Obj(oxName=self.oxidizer, fuelName=self.fuel)
        
        # State variables for network compatibility
        self.P = Pa 
        self.T = 293.15
        self.m = 0.0
        self.V = 0.001
        self.d = 1.2
        self.Q = None
        self.fill_level = 0.0
        self.fluid = "CombustionGas"
        
        # Engine performance states
        self.mdot_ox = 0.0
        self.mdot_fu = 0.0
        self.MR = 0.0
        self.cstar = 0.0
        self.thrust = 0.0
        self.Isp = 0.0
        
        # Expanded history dict
        self.history = {k: [] for k in [
            "time", "P", "T", "m", "d", "Q", "fill_level", 
            "mdot_ox", "mdot_fu", "MR", "cstar", "thrust", "Isp"
        ]}

    def update(self, mdot_net, Hdot_net, dt):
        self.mdot_ox = abs(self.ox_conn.mdot)
        self.mdot_fu = abs(self.fuel_conn.mdot)
        mdot_total = self.mdot_ox + self.mdot_fu
        
        if mdot_total < 0.01:
            self.P, self.T, self.MR, self.thrust, self.Isp, self.cstar = self.Pa, 293.15, 0.0, 0.0, 0.0, 0.0
            return
            
        self.P = self.P
        
        if self.mdot_fu < 1e-4: self.MR = 99.0 
        elif self.mdot_ox < 1e-4: self.MR = 0.01
        else: self.MR = self.mdot_ox / self.mdot_fu
        
        MR_cea = max(0.5, min(self.MR, 10.0))
        Pc_psia = max(self.P / 6894.75729, 14.7)
        
        # Evaluate dynamic efficiencies
        curr_eta_cstar = self.eta_cstar(Pc_psia) if callable(self.eta_cstar) else self.eta_cstar
        curr_eta_cf = self.eta_cf(Pc_psia) if callable(self.eta_cf) else self.eta_cf
        
        try:
            cstar_ft = self.cea.get_Cstar(Pc=Pc_psia, MR=MR_cea)
            cstar_ideal = cstar_ft * 0.3048 
            self.T = self.cea.get_Tcomb(Pc=Pc_psia, MR=MR_cea) * (5.0 / 9.0)
        except Exception:
            cstar_ideal, self.T = 1000.0, 1000.0

        self.cstar = cstar_ideal * curr_eta_cstar
        Pc_target = (mdot_total * self.cstar) / self.At
        
        relaxation = 0.3 
        self.P = max((self.P * (1.0 - relaxation)) + (Pc_target * relaxation), self.Pa)

        if self.P > (self.Pa * 1.5):
            eps = self.Ae / self.At
            Pc_psia = self.P / 6894.75729
            try:
                isp_amb = self.cea.estimate_Ambient_Isp(Pc=Pc_psia, MR=MR_cea, eps=eps, Pamb=(self.Pa/6894.75729))[0]
                # Total Isp efficiency is roughly eta_cstar * eta_cf
                self.Isp = isp_amb * curr_eta_cstar * curr_eta_cf 
                self.thrust = self.Isp * 9.81 * mdot_total
            except Exception:
                self.thrust, self.Isp = 0.0, 0.0
        else:
            self.thrust, self.Isp = 0.0, 0.0

        self.P = self.P

    def log_state(self, t=0.0):
        self.history["time"].append(t)
        self.history["P"].append(self.P)
        self.history["T"].append(self.T)
        self.history["m"].append(self.m) # Dummy for plotter
        self.history["d"].append(self.d) # Dummy for plotter
        self.history["Q"].append(self.Q) # Dummy for plotter
        self.history["fill_level"].append(self.fill_level) # Dummy for plotter
        
        self.history["mdot_ox"].append(self.mdot_ox)
        self.history["mdot_fu"].append(self.mdot_fu)
        self.history["MR"].append(self.MR)
        self.history["cstar"].append(self.cstar)
        self.history["thrust"].append(self.thrust)
        self.history["Isp"].append(self.Isp)


class Tank(Node):
    """
    Two-phase Tank Node with Heat Transfer.
    The 'Tank' instance itself represents the Liquid node (bottom).
    It contains a .ullage attribute which is the Gas node (top).
    
    The two nodes are coupled by a Volume constraint: V_liq + V_gas = V_tank.
    Pressure is iterated until this constraint is met.
    
    Includes lumped capacitance heat transfer between Liquid and Ullage (Collapse).
    """
    def __init__(self, V_total_L, fluid_liq, m_liq, T_liq, fluid_ullage, P_ullage, T_ullage, 
                 radius=None, htc=50.0, name="tank"):
        """
        Args:
            V_total_L: Total tank volume in Liters.
            fluid_liq: String name of liquid fluid.
            m_liq: Initial mass of liquid (kg).
            T_liq: Initial temperature of liquid (K).
            fluid_ullage: String name of ullage fluid.
            P_ullage: Initial ullage pressure (Pa).
            T_ullage: Initial ullage temperature (K).
            radius (optional): Tank radius (m) for heat transfer area calc. 
            htc (optional): Heat transfer coefficient (W/m^2K) for ullage-liquid HT.
        """
        self.V_total = float(V_total_L) / 1000.0  # Store fixed tank volume [m^3]
        self.radius = radius
        self.htc = htc
        
        # Calculate Heat Transfer Area
        if self.radius:
            self.area_interface = np.pi * self.radius**2
        else:
            self.area_interface = 0.0

        # --- 1. Calculate Initial Conditions to Satisfy Inputs ---
        
        # Calculate Liquid Density/Volume based on Ullage Pressure (assuming mechanical equilibrium)
        # Note: Neglecting hydrostatic head for 0D initialization
        try:
            rho_liq = PropsSI_auto('D', 'P', P_ullage, 'T', T_liq, fluid_liq)
        except:
            # Fallback if P/T combo is invalid (e.g. subcooled logic fail), try sat liquid
            rho_liq = PropsSI_auto('D', 'P', P_ullage, 'Q', 0, fluid_liq)

        V_liq = m_liq / rho_liq
        
        # Calculate Ullage Volume and Mass
        V_ull = self.V_total - V_liq
        if V_ull < 0:
            raise ValueError(f"Tank {name}: Liquid mass {m_liq}kg exceeds total volume at {T_liq}K.")
            
        rho_ull = PropsSI_auto('D', 'P', P_ullage, 'T', T_ullage, fluid_ullage)
        m_ull = V_ull * rho_ull

        # --- 2. Initialize Nodes ---
        
        # Initialize Liquid Node (Self)
        # We assume the liquid fills its computed volume initially
        super().__init__(fluid_liq, m_liq, V_liq * 1000.0, T_liq, name=name)
        
        # Initialize Ullage Node
        self.ullage = Node(fluid_ullage, m_ull, V_ull * 1000.0, T_ullage, name=f"{name}_ullage")
        
        # Correct Fill Level for Tank Geometry (Node class defaults fill_level to phase fraction)
        self.fill_level = self.V / self.V_total

        # --- 3. Force initial balance ---
        self._balance_volumes()

    def update(self, mdot_l, Hdot_l, mdot_g, Hdot_g, dt):
        """
        Custom update that handles mass/energy fluxes for both phases,
        applies heat transfer, and enforces the shared pressure/volume constraint.
        """
        # --- Step 0: Heat Transfer (Ullage Collapse / Evap) ---
        # Lumped capacitance: Q = h * A * (T_ullage - T_liquid)
        if self.area_interface > 0:
            # Positive Q flows from Ullage -> Liquid
            Q_transfer = self.htc * self.area_interface * (self.ullage.T - self.T)
            
            # Hdot_l += Q_transfer
            Hdot_g -= Q_transfer

        # --- Step 1: Integrate Mass and Energy (Euler Step) ---
        # Expansion Work Transfer
        # The Ullage expands and does work on the Liquid: W = P * dV
        rho_l = self.d_l if (hasattr(self, 'd_l') and self.d_l > 1) else self.d
        if rho_l < 1: rho_l = 1000.0 # Safety div/0

        # Calculate Work Rate (Watts) = P * (dV_ullage/dt)
        # If liquid drains (mdot_l < 0), Ullage expands (dV_ullage > 0).
        # Work done BY ullage is positive.
        vol_change_rate = -(mdot_l / rho_l) 
        work_power = self.P * vol_change_rate

        # Apply Energy Transfer
        # Ullage DOES work (Loses Energy)
        Hdot_g -= work_power
        # Liquid RECEIVES work (Gains Energy - which is then exported as Flow Work Pv)
        Hdot_l += work_power

        self.m += mdot_l * dt
        self.U += Hdot_l * dt
        self.ullage.m += mdot_g * dt
        self.ullage.U += Hdot_g * dt
        # Numerical safety floors
        if self.m < 1e-12: self.m = 1e-12
        if self.ullage.m < 1e-12: self.ullage.m = 1e-12

        # --- Step 2: Solve for shared Pressure ---
        # This function adjusts self.V and self.ullage.V until they sum to V_total
        self._balance_volumes()

        # --- Step 3: Update Fluid States (Flash) ---
        # Update Liquid State
        self.d = self.m / self.V
        self._flash_from_DU(self.d, self.U)
        
        # Override fill_level to mean "Tank Fill" (not liquid phase fraction)
        self.fill_level = self.V / self.V_total
        # Update Ullage State
        self.ullage.d = self.ullage.m / self.ullage.V
        self.ullage._flash_from_DU(self.ullage.d, self.ullage.U)

    def _balance_volumes(self):
        """
        Iteratively find the Pressure P such that:
        Volume_Liq(P, u_l) + Volume_Gas(P, u_g) = V_Total
        Uses Internal Energy (u) because U is conserved, whereas H varies with P.
        """
        # Specific internal energies (held constant during P solve)
        u_l = self.U / self.m
        u_g = self.ullage.U / self.ullage.m

        # Initial Guess for P
        p_guess = self.P
        p_step = 1000.0 # 1 kPa perturbation
        
        # Pressure balance solver loop (Secant Method)
        for i in range(20):
            # Calculate Residual for P1
            err1 = self._get_vol_error(p_guess, u_l, u_g)
            
            if abs(err1) < 1e-6: # Tolerance: 1 mL
                break
                
            # Calculate Residual for P2 (perturbed)
            p_guess_2 = p_guess + p_step
            err2 = self._get_vol_error(p_guess_2, u_l, u_g)
            
            # Secant update
            denom = (err2 - err1)
            if abs(denom) < 1e-12:
                break # Jacobian singular, stick with current P
            
            p_new = p_guess - err1 * (p_step / denom)
            
            # Safety clamp for pressure (non-negative, min 100 Pa)
            if p_new < 100: p_new = 100.0 
            
            p_guess = p_new
            
        # Apply the final Volumes based on the found P
        # NOTE: Changed 'E' to 'U' for CoolProp/Standard compatibility
        rho_l = PropsSI_auto('D', 'P', p_guess, 'E', u_l, self.fluid)
        rho_g = PropsSI_auto('D', 'P', p_guess, 'E', u_g, self.ullage.fluid)
        
        self.V = self.m / rho_l
        self.ullage.V = self.ullage.m / rho_g
        self.P = p_guess
        self.ullage.P = p_guess

    def _get_vol_error(self, p, u_l, u_g):
        """ 
        Helper to calculate Volume Error at a given Pressure.
        Updated to take (u) instead of (h).
        """
        try:
            # Get densities at candidate Pressure & Fixed Internal Energy
            # NOTE: Changed 'E' to 'U'
            rho_l = PropsSI_auto('D', 'P', p, 'E', u_l, self.fluid)
            rho_g = PropsSI_auto('D', 'P', p, 'E', u_g, self.ullage.fluid)
            
            v_l = self.m / rho_l
            v_g = self.ullage.m / rho_g
            
            return (v_l + v_g) - self.V_total
        except:
            # If PropsSI fails (e.g. out of bounds), return large error
            return 1.0
        

# ==============================================================================
# CONNECTION CLASS AND SUBCLASSES
# ==============================================================================
"""
TODO: 
REDEFINE CONNECTION CLASSES TO MODEL INERTANCE
ADD INLET AND OUTLET LOCATIONS FOR CONNECTION (IDK IF THIS ACTUALLY MATTERS)
REDEFINE QDOT IN WITH HEATTRANSFER OBJECT
REDEFINE ALL CONTROL USING CONTROL OBJECT/SUBCLASSES
"""

class Connection():
    """
    Connection class. Defined by CdA (m^2), qdot (J/s), location on node (0-1), and state (open, closed).
    Initialized by CdA, qdot, location, and normal state.
    """
    def __init__(self, CdA, qdot=0.0, location=0.0, normal_state=True, checking=True, name="connection"):
        self.CdA = CdA
        self.name = name
        self.dP = 0 # to be used
        self.qdot = qdot
        self.location = location  # normalized height 0-1
        self.state = normal_state
        self.checking = checking
        self.mdot = 0
        self.Hdot = 0
        self.Q = None
        self.history = {k: [] for k in ["time","CdA", "qdot","state", "mdot","Hdot","dP", "Q"]}

    def mdot_Hdot(self, node1, node2):
        """
        Return mdot (kg/s), Hdot (J/s) where positive means mass/enthalpy flows node1 -> node2.
        Includes Dyer model for two-phase flow (flashing).
        """
        # Check if connection is open
        if not self.state:
            return 0.0, 0.0

        dP = node1.P - node2.P
        if self.checking and dP < 0:
            return 0.0, 0.0
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        # Determine donor and receiver
        if dP > 0:
            donor, receiver = node1, node2
        else:
            donor, receiver = node2, node1

        # --- Phase / Source Selection Logic ---
        # Check if donor is a Tank (has ullage)
        if hasattr(donor, 'ullage'):
            # It is a Tank: Compare connection location to tank fill level
            if self.location > donor.fill_level:
                # Connection is in the Ullage (Gas)
                source_node = donor.ullage
            else:
                # Connection is in the Liquid
                source_node = donor
        else:
            # Standard Node: Logic relies on internal phase fraction if two-phase
            if donor.fill_level > self.location:
                source_node = donor # Liquid part of node (if split) or Bulk
            else:
                source_node = donor # Gas part... 
                # Note: For standard nodes, we often just use bulk properties unless we specifically 
                # implemented separated h_l/h_v access. Below we handle the bulk/separation.

        # Retrieve source properties
        # If we selected the Ullage node, it acts as a single-phase gas node.
        # If we selected the Liquid node (Tank), it acts as a liquid node (potentially 2-phase if boiling).
        
        # Are we pulling liquid or vapor?
        # If source_node is ullage -> It's gas.
        # If source_node is tank/node -> Check its internal phase.
        
        # Simplified Logic using available properties on the chosen source node:
        # If the source node is "pure" (like the ullage node), h_l approx h_v approx h.
        
        if source_node.fill_level > self.location: # Mostly liquid
            h_stream = source_node.h_l if hasattr(source_node, 'h_l') else source_node.h
            d_stream = source_node.d_l if hasattr(source_node, 'd_l') else source_node.d
        else:
            h_stream = source_node.h_v if hasattr(source_node, 'h_v') else source_node.h
            d_stream = source_node.d_v if hasattr(source_node, 'd_v') else source_node.d

        # For explicit Tank Ullage access, we override the above:
        if hasattr(donor, 'ullage') and self.location > donor.fill_level:
            # We are explicitly in the ullage node
            h_stream = donor.ullage.h
            d_stream = donor.ullage.d
        elif hasattr(donor, 'ullage') and self.location <= donor.fill_level:
            # We are explicitly in the liquid node
            h_stream = donor.h
            d_stream = donor.d
        abs_dP = abs(dP)
        self.dP = abs(dP)  # logging
        # Determine phase for flow model
        donor_phase = CP.PhaseSI('D', source_node.d, 'H', source_node.h, source_node.fluid)
        # --- GAS/CHOKED FLOW --- 
        if donor_phase in ("gas", "supercritical", "supercritical_gas"):
            gamma = source_node.gamma
            R = source_node.R
            Tdon = source_node.T
            crit_factor = ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
            Pcrit = source_node.P * crit_factor

            if receiver.P > Pcrit:
                # Unchoked subsonic gas flow
                mdot_mag = self.CdA * source_node.P * np.sqrt(2 * abs(1 - (receiver.P / source_node.P) ** ((gamma - 1) / gamma)) / (R * Tdon))
            else:
                # Choked
                mdot_mag = self.CdA * source_node.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
            Hdot = mdot_mag * h_stream
        
        # --- TWO-PHASE (Dyer model) ---
        elif donor_phase == "twophase":
            h_liq = PropsSI_auto('H', 'P', receiver.P, 'Q', 0, source_node.fluid)
            h_vap = PropsSI_auto('H', 'P', receiver.P, 'Q', 1, source_node.fluid)
            Pv = PropsSI_auto('P', 'T', source_node.T, 'Q', 1, source_node.fluid)

            # Single-phase incompressible term (SPI)
            mdot_spi = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * abs_dP)

            # Homogeneous equilibrium model term (HEM)
            try:
                h1 = h_stream
                h2 = PropsSI_auto('H', 'P', receiver.P, 'S', source_node.s, source_node.fluid)
                rho2p = PropsSI_auto('D', 'P', receiver.P, 'S', source_node.s, source_node.fluid)
                mdot_hem = self.CdA * rho2p * np.sqrt(2.0 * max(h1 - h2, 1e-9))
            except Exception:
                mdot_hem = mdot_spi

            # Dyer blending factor
            r = 1  # tunable, change based on test data
            kappa = r * (source_node.P - receiver.P) / max(Pv - receiver.P, 1e-6) # can also manually set kappa (2 is a good conservative value)

            # Dyer blended mass flow
            mdot_mag = (kappa / (1 + kappa)) * mdot_spi + (1 / (1 + kappa)) * mdot_hem

            # Enthalpy flow rate
            self.Q = PropsSI_auto('Q', 'P', receiver.P, 'H', h_stream, source_node.fluid)
            if 0 <= self.Q <= 1:
                Hdot = mdot_mag * (self.Q * h_vap + (1 - self.Q) * h_liq)
            else:
                Hdot = mdot_mag * h_stream

        # --- LIQUID ---
        else:
            mdot_mag = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * abs_dP)
            Hdot = mdot_mag * h_stream
        
        # Sign convention
        if donor is node1:
            mdot = mdot_mag
        else:
            mdot = -mdot_mag
            Hdot = - Hdot

        Hdot += self.qdot  # add any heat leak term
        self.mdot, self.Hdot = mdot, Hdot
        return mdot, Hdot

    def log_state(self, t=0.0):
        """
        Log node state at each timestep throughout a network sim.
        """
        self.history["time"].append(t)
        self.history["CdA"].append(self.CdA)
        self.history["qdot"].append(self.qdot)
        self.history["state"].append(self.state)
        self.history["mdot"].append(self.mdot)
        self.history["Hdot"].append(self.Hdot)
        self.history["dP"].append(self.dP)
        self.history["Q"].append(self.Q)


class Line(Connection):
    """
    Line class (subclass of Connection). Defined by physical pipe geometry.
    Dynamically updates CdA based on Darcy-Weisbach (liquids) or True Fanno Flow (gases).
    """
    def __init__(self, ID, length, roughness, qdot=0.0, location=0.0, normal_state=True, checking=True, name="line"):
        self.ID = ID
        self.length = length
        self.roughness = roughness
        self.Area = np.pi * (self.ID**2) / 4.0
        
        if self.roughness > 0:
            f_init = (-2.0 * np.log10(self.roughness / (3.7 * self.ID)))**-2
        else:
            f_init = 0.015 
        
        #dp = 0.5 * (self.length / self.ID) * 
        
        initial_CdA = self.Area / np.sqrt(f_init * self.length / self.ID)
        
        super().__init__(initial_CdA, qdot, location, normal_state, checking, name)

    def mdot_Hdot(self, node1, node2):
        if not self.state:
            return 0.0, 0.0

        dP = node1.P - node2.P
        if self.checking and dP < 0:
            return 0.0, 0.0
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        donor = node1 if dP > 0 else node2
        receiver = node2 if dP > 0 else node1
        
        if hasattr(donor, 'ullage') and self.location > donor.fill_level:
            source_node = donor.ullage
        else:
            source_node = donor
            
        if source_node.fill_level > self.location:
            h_stream = source_node.h_l if hasattr(source_node, 'h_l') else source_node.h
            d_stream = source_node.d_l if hasattr(source_node, 'd_l') else source_node.d
        else:
            h_stream = source_node.h_v if hasattr(source_node, 'h_v') else source_node.h
            d_stream = source_node.d_v if hasattr(source_node, 'd_v') else source_node.d

        try:
            mu = PropsSI_auto('V', 'D', d_stream, 'H', h_stream, source_node.fluid)
        except Exception:
            mu = 1e-4 

        if abs(self.mdot) < 1e-12:
            V_guess = np.sqrt(2.0 * abs(dP) / max(d_stream, 1e-6))
            mdot_guess = d_stream * self.Area * V_guess
        else:
            mdot_guess = abs(self.mdot)

        Re = (4.0 * mdot_guess) / (np.pi * self.ID * mu)


        # Smooth Transitional Blending to prevent Numerical Chatter
        if Re < 2000:
            f = 64.0 / max(Re, 1e-6)
        elif Re > 3000:
            f = (-1.8 * np.log10((self.roughness / (3.7 * self.ID))**1.11 + 6.9 / Re))**-2
        else:
            # Cosine interpolation between Laminar and Turbulent
            f_lam = 64.0 / 2000.0
            f_turb = (-1.8 * np.log10((self.roughness / (3.7 * self.ID))**1.11 + 6.9 / 3000.0))**-2
            weight = (1.0 - np.cos(np.pi * (Re - 2000.0) / 1000.0)) / 2.0
            f = f_lam * (1.0 - weight) + f_turb * weight

        donor_phase = CP.PhaseSI('D', d_stream, 'H', h_stream, source_node.fluid)
        
        # --- TRUE FANNO FLOW LOGIC ---
        if donor_phase in ("gas", "supercritical", "supercritical_gas"):
            gamma = source_node.gamma
            R = source_node.R
            fLD = f * self.length / self.ID
            
            # 1. Fanno Parameter Function F(M)
            def F_param(M):
                if M <= 1e-5: return 1e9
                if M >= 1.0: return 0.0
                term1 = (1.0 - M**2) / (gamma * M**2)
                term2 = ((gamma + 1.0) / (2.0 * gamma)) * np.log(((gamma + 1.0) * M**2) / (2.0 + (gamma - 1.0) * M**2))
                return term1 + term2
                
            # 2. Find Choked Inlet Mach Number (M1) using Bisection
            low, high = 1e-5, 1.0
            for _ in range(25): # 25 iterations is highly accurate and ultra-fast
                mid = (low + high) / 2.0
                if F_param(mid) > fLD:
                    low = mid
                else:
                    high = mid
            M1_choked = (low + high) / 2.0
            
            # 3. Calculate Static Inlet Properties (Assuming Node P, T are Stagnation)
            T1_static = source_node.T / (1.0 + 0.5 * (gamma - 1.0) * M1_choked**2)
            P1_static = source_node.P / (1.0 + 0.5 * (gamma - 1.0) * M1_choked**2)**(gamma / (gamma - 1.0))
            
            # 4. Calculate Exit Pressure if Choked (M2 = 1)
            P2_choked = P1_static * M1_choked * np.sqrt((2.0 + (gamma - 1.0) * M1_choked**2) / (gamma + 1.0))
            
            # 5. Determine Choked vs Unchoked and Assign CdA
            if receiver.P <= P2_choked:
                # CHOKED: Calculate true Fanno mass flow
                V1 = M1_choked * np.sqrt(gamma * R * T1_static)
                rho1 = P1_static / (R * T1_static)
                mdot_fanno = rho1 * self.Area * V1
                
                # Back-calculate CdA to force the base class to yield mdot_fanno
                crit_factor = ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
                base_mdot_1CdA = (source_node.P / np.sqrt(max(source_node.T, 1e-8))) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
                self.CdA = mdot_fanno / base_mdot_1CdA
            else:
                # UNCHOKED: Use the compressible equivalent Area mapped from fLD
                # (Iterative unchoked Fanno solving creates instability in transient networks, 
                # this equivalent resistance method is standard for 1D unchoked solvers)
                self.CdA = self.Area / np.sqrt(1.0 + fLD)
                
        # --- INCOMPRESSIBLE / TWO-PHASE LOGIC ---
        else:
            self.CdA = self.Area / np.sqrt(f * self.length / self.ID)

        # Delegate the actual mass/enthalpy flow calculation to the parent class
        return super().mdot_Hdot(node1, node2)
    

class Series(Connection):
    """
    Class to combine all the connections in series between nodes.
    Features an active flow solver for components in "mdot" mode.
    """
    def __init__(self, connections, name='series'):
        self.connections = connections
        normal_state = all(c.state for c in self.connections)
        series_checking = any(getattr(c, 'checking', False) for c in self.connections)
        super().__init__(0.0, 0.0, connections[0].location, normal_state, series_checking, name)
        
        # Trigger an initial calculation of CdA
        self._recalc_cda()

    def _recalc_cda(self):
        sum_inverse_squares = 0.0
        for c in self.connections:
            eff_state = getattr(c, 'state', True)
            eff_CdA = c.CdA if eff_state else 0.0
            if eff_CdA > 0:
                sum_inverse_squares += 1.0 / (eff_CdA**2)
            else:
                sum_inverse_squares = float('inf')
                break
                
        if sum_inverse_squares == float('inf') or sum_inverse_squares == 0:
            self.CdA = 0.0
        else:
            self.CdA = 1.0 / np.sqrt(sum_inverse_squares)

    def mdot_Hdot(self, node1, node2):
        total_dP = abs(node1.P - node2.P)
        
        if getattr(self, 'state', True) == False:
            self.CdA, self.mdot, self.Hdot = 0.0, 0.0, 0.0
            for c in self.connections:
                c.mdot, c.Hdot = 0.0, 0.0
                c.dP = total_dP if getattr(c, 'state', True) == False else 0.0
            return 0.0, 0.0

        # 1. Update passive components and identify Active Controllers
        mfc_components = []
        for c in self.connections:
            if getattr(c, 'mode', '') == 'mdot':
                c.state = max(0.0, min(1.0, float(c.state)))
                c.target_mdot = c.CdA_max * c.state
                mfc_components.append(c)
            else:
                c.mdot_Hdot(node1, node2) # Lines update f, Re here

        # 2. Active Flow Solver
        if mfc_components:
            mfc = mfc_components[0] # Assuming one active valve per series
            if mfc.target_mdot == 0.0:
                mfc.CdA = 0.0
            else:
                # Calculate resistance of all OTHER components
                sum_inv_sq_other = 0.0
                for c in self.connections:
                    if c is not mfc:
                        eff_state = getattr(c, 'state', True)
                        eff_CdA = c.CdA if eff_state else 0.0
                        if eff_CdA > 0:
                            sum_inv_sq_other += 1.0 / (eff_CdA**2)
                        else:
                            sum_inv_sq_other = float('inf')
                            break
                
                if sum_inv_sq_other == float('inf'):
                    mfc.CdA = 0.0
                else:
                    # Find required CdA for the WHOLE series
                    temp_CdA = self.CdA
                    self.CdA = 1.0
                    flux_mdot, _ = super().mdot_Hdot(node1, node2)
                    self.CdA = temp_CdA # restore
                    
                    if flux_mdot == 0:
                        mfc.CdA = 0.0
                    else:
                        cda_req_series = mfc.target_mdot / abs(flux_mdot)
                        inv_sq_req = 1.0 / (cda_req_series**2)
                        
                        # Apply Inverse Square Resistance Law
                        inv_sq_valve = inv_sq_req - sum_inv_sq_other
                        
                        if inv_sq_valve <= 0:
                            # Other components are too restrictive. Fully open valve.
                            mfc.CdA = mfc.CdA_max
                        else:
                            mfc.CdA = 1.0 / np.sqrt(inv_sq_valve)
                            if mfc.CdA > mfc.CdA_max:
                                mfc.CdA = mfc.CdA_max

        # 3. Calculate final equivalent CdA dynamically
        self._recalc_cda()
        
        if self.CdA == 0.0:
            self.mdot, self.Hdot = 0.0, 0.0
            for c in self.connections:
                c.mdot, c.Hdot = 0.0, 0.0
                eff_state = getattr(c, 'state', True)
                c.dP = total_dP if (c.CdA == 0 or not eff_state) else 0.0
            return 0.0, 0.0

        # 4. Total Flow
        mdot, Hdot = super().mdot_Hdot(node1, node2)

        # 5. Propagate back down
        for c in self.connections:
            c.mdot = mdot
            c.Hdot = Hdot
            eff_state = getattr(c, 'state', True)
            eff_CdA = c.CdA if eff_state else 0.0
            
            if self.CdA > 0 and eff_CdA > 0:
                c.dP = total_dP * ((self.CdA / eff_CdA)**2)
            elif eff_CdA == 0.0:
                c.dP = total_dP
            else:
                c.dP = 0.0

        return mdot, Hdot

    def log_state(self, t=0.0):
        super().log_state(t)
        for c in self.connections:
            c.log_state(t)

class Regulator(Connection):
    def __init__(self, CdA, set_pressure, droop_curve=None, qdot=0.0, location=0.0, normal_state=True):
        """
        NOTE: STILL IN DEVELOPLENT
        A pressure regulator connection that limits downstream pressure. Defined by: CdA (m^2), set_pressure (Pa),
        droop_curve (function that maps mdot -> pressure drop (Pa)), qdot (J/s), location (0-1), state (open, closed).
        """
        super().__init__(CdA, qdot, location, normal_state)
        self.set_pressure = set_pressure
        self.droop_curve = droop_curve  # function handle: f(mdot) -> ΔP droop

    def mdot_Hdot(self, node1, node2):
        """
        Computes mdot and Hdot across the regulator.
        The regulator limits downstream pressure to set_pressure (minus droop if defined).
        Args:
            node1, node1 (Node): nodes connected by this connection
        """
        if not self.state:
            return 0.0, 0.0

        # Determine upstream and downstream
        dP = node1.P - node2.P
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        if dP > 0:
            upstream, downstream = node1, node2
        else:
            upstream, downstream = node2, node1

        # Target downstream pressure
        P_down_target = self.set_pressure

        # Apply droop curve if defined
        if self.droop_curve is not None:
            # iterative droop correction: assume mdot ≈ previous mdot, or start with 0
            # droop curve returns positive ΔP loss at higher flows
            P_down_target -= self.droop_curve(abs(dP))  

        # Clamp downstream pressure to not exceed target
        if downstream.P < P_down_target:
            # regulator closed: no flow (receiver pressure too low)
            return 0.0, 0.0
        else:
            # regulator open: limit flow so that downstream ≈ setpoint
            effective_dP = max(upstream.P - P_down_target, 0.0)

        # Now use inherited orifice logic for the flow
        # Source Selection Logic (Regulator specific)
        if hasattr(upstream, 'ullage') and self.location > upstream.fill_level:
             source = upstream.ullage
        else:
             source = upstream

        if source.fill_level > self.location:
            h_stream = source.h_l if hasattr(source, 'h_l') else source.h
            d_stream = source.d_l if hasattr(source, 'd_l') else source.d
        else:
            h_stream = source.h_v if hasattr(source, 'h_v') else source.h
            d_stream = source.d_v if hasattr(source, 'd_v') else source.d
            
        # Simplified access
        h_stream = source.h
        d_stream = source.d

        donor_phase = CP.PhaseSI('D', source.d, 'H', source.h, source.fluid)
        if donor_phase in ("gas", "supercritical") and source.Cp and source.Cv and source.R:
            gamma = source.gamma
            R = source.R
            Tdon = source.T
            crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
            mdot_mag = self.CdA * source.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
        else:
            mdot_mag = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * effective_dP)

        # Sign convention: positive mdot -> node1 -> node2
        mdot = mdot_mag if upstream is node1 else -mdot_mag

        Hdot = mdot * h_stream + self.qdot
        return mdot, Hdot


class BangBang(Connection):
    """
    Bang-bang valve controller.
    Opens if downstream P < target - hysteresis.
    Closes if downstream P > target + hysteresis.
    """
    def __init__(self, CdA, target_pressure, hysteresis=0.0, qdot=0.0, location=0.0, normal_state=True, checking=True, name="bang_bang"):
        super().__init__(CdA, qdot, location, normal_state, checking, name)
        self.target_pressure = target_pressure
        self.hysteresis = hysteresis

    def update_control(self, node1, node2):
        """
        Determines the downstream node and updates state based on pressure.
        """
        # Determine which node is downstream based on pressure gradient
        # (Or you could enforce a direction, but this is more general)
        if node1.P > node2.P:
            downstream = node2
        else:
            downstream = node1

        # Control Logic with Hysteresis
        if self.state:
            # If currently OPEN, stay open until we hit the upper limit
            if downstream.P > (self.target_pressure + self.hysteresis):
                self.state = False
        else:
            # If currently CLOSED, stay closed until we hit the lower limit
            if downstream.P < (self.target_pressure - self.hysteresis):
                self.state = True


class SharpEdgedOrifice(Connection):
    """
    Subclass of Node to represent a sharp-edged orifice.
    """
    # TODO
    pass


class ThrottleValve(Connection):
    """
    Proportional throttle valve.
    mode="cda":  'state' (0-1) linearly scales the flow area up to CdA_max (m^2).
    mode="mdot": 'state' (0-1) linearly scales the mass flow up to CdA_max (which acts as Max Target mdot in kg/s).
    """
    def __init__(self, CdA_max, qdot=0.0, location=0.0, normal_state=0.0, checking=True, mode="cda", name="throttle_valve"):
        self.CdA_max = CdA_max
        self.mode = mode.lower()
        self.target_mdot = 0.0
        initial_state = max(0.0, min(1.0, float(normal_state)))
        initial_CdA = self.CdA_max * initial_state if self.mode == "cda" else 0.0
        
        super().__init__(initial_CdA, qdot, location, initial_state, checking, name)

    def mdot_Hdot(self, node1, node2):
        self.state = max(0.0, min(1.0, float(self.state)))
        
        if self.state == 0.0:
            self.mdot, self.Hdot, self.CdA = 0.0, 0.0, 0.0
            return 0.0, 0.0

        if self.mode == "cda":
            self.CdA = self.CdA_max * self.state
            return super().mdot_Hdot(node1, node2)
            
        elif self.mode == "mdot":
            self.target_mdot = self.CdA_max * self.state
            self.CdA = 1.0
            mdot_flux, Hdot_flux = super().mdot_Hdot(node1, node2)
            
            if mdot_flux == 0.0:
                self.CdA, self.mdot, self.Hdot = 0.0, 0.0, 0.0
                return 0.0, 0.0
                
            self.CdA = self.target_mdot / abs(mdot_flux)
            
            # Clamp to physical limits (Cannot open wider than 100%)
            if self.CdA > self.CdA_max:
                self.CdA = self.CdA_max
                
            self.mdot = mdot_flux * self.CdA
            Hdot_fluid = Hdot_flux - self.qdot
            self.Hdot = (Hdot_fluid * self.CdA) + self.qdot
            
            return self.mdot, self.Hdot
                        
    
# ==============================================================================
# NETWORK CLASS AND PLOTTING
# ==============================================================================

class Network():
    """
    Network class. Defined by a graph of connections and nodes.
    Automatically detects 'Tank' objects to handle coupled liquid/ullage updates.
    """
    def __init__(self, graph):
        self.graph = graph  # {connection: (node1, node2)}
        
        # Pre-scan the graph to identify Tank objects.
        self.tanks = set()
        for pair in self.graph.values():
            for node in pair:
                if type(node).__name__ == 'Tank': 
                    self.tanks.add(node)

    def sim(self, t, dt, actions={}, verbose_steps=5):
        """
        Runs transient sim. Handles standard Nodes, coupled Tank Nodes, and Active Valves.
        Includes Smart Routing to direct flux to/from the correct Tank phase (Liquid vs Ullage).
        """
        steps = int(t / dt)
        
        for i in range(steps):
            time_now = round(i * dt, 1)
            
            # 1. Apply Scripted Actions
            if time_now in actions:
                # Iterate through all scheduled events for this specific timestep
                for conn, state in actions[time_now]:
                    conn.state = state
                    if verbose_steps > 0:
                        print(f"--- Action at {time_now}s: {conn.name} set to {state} ---")

            # 2. Update Active Components (BangBang, Regulators)
            for conn, (n1, n2) in self.graph.items():
                # Check top-level connections
                if hasattr(conn, 'update_control'):
                    conn.update_control(n1, n2)
                
                # Check encapsulated connections (like inside a Series object)
                if hasattr(conn, 'connections'):
                    for sub_conn in conn.connections:
                        if hasattr(sub_conn, 'update_control'):
                            # Pass the global nodes so the smart valve can read the global dP
                            sub_conn.update_control(n1, n2)

            # 3. Initialize Flux Containers
            # Must include ALL base nodes AND their sub-nodes (ullage)
            all_nodes = set()
            for pair in self.graph.values():
                for node in pair:
                    all_nodes.add(node)
                    if hasattr(node, 'ullage'):
                        all_nodes.add(node.ullage)
            
            mdot_contrib = {n: 0.0 for n in all_nodes}
            Hdot_contrib = {n: 0.0 for n in all_nodes}

            # 4. Compute and Route Fluxes
            for conn, (n1, n2) in self.graph.items():
                mdot, Hdot = conn.mdot_Hdot(n1, n2)
                
                # --- SMART ROUTING LOGIC ---
                # Determine the effective Source Node
                if hasattr(n1, 'ullage') and hasattr(n1, 'fill_level'):
                    # Node 1 is a Tank: Route based on connection location
                    effective_n1 = n1.ullage if conn.location > n1.fill_level else n1
                else:
                    effective_n1 = n1

                # Determine the effective Target Node
                if hasattr(n2, 'ullage') and hasattr(n2, 'fill_level'):
                    # Node 2 is a Tank: Route based on connection location
                    effective_n2 = n2.ullage if conn.location > n2.fill_level else n2
                else:
                    effective_n2 = n2

                # Apply Fluxes to the EFFECTIVE nodes
                # Flow convention: n1 -> n2 is positive
                mdot_contrib[effective_n1] -= mdot
                mdot_contrib[effective_n2] += mdot
                Hdot_contrib[effective_n1] -= Hdot
                Hdot_contrib[effective_n2] += Hdot
                
                conn.log_state(time_now)

            # 5. Update Nodes
            processed_nodes = set()

            # --- A. Update Tanks (Coupled Liquid + Ullage) ---
            for tank in self.tanks:
                # Tank.update expects fluxes for both liquid and ullage separately
                mdot_l = mdot_contrib.get(tank, 0.0)
                Hdot_l = Hdot_contrib.get(tank, 0.0)
                mdot_g = mdot_contrib.get(tank.ullage, 0.0)
                Hdot_g = Hdot_contrib.get(tank.ullage, 0.0)

                tank.update(mdot_l, Hdot_l, mdot_g, Hdot_g, dt)
                processed_nodes.add(tank)
                processed_nodes.add(tank.ullage)
                
                if i < verbose_steps:
                    print(f"[t={time_now:.4f}] {tank.name} P={tank.P/1e5:.2f} bar, Fill={tank.fill_level:.2f}")

            # --- B. Update Standard Nodes ---
            for node in mdot_contrib:
                if node not in processed_nodes:
                    # Only update if there was flux (optimization) or if it's an active node
                    if abs(mdot_contrib[node]) > 0 or abs(Hdot_contrib[node]) > 0 or not hasattr(node, 'update'):
                         # Note: Ambient nodes have empty update() pass, so it's safe
                        node.update(mdot_contrib[node], Hdot_contrib[node], dt)
                    processed_nodes.add(node)
                    
                    if i < verbose_steps and abs(mdot_contrib[node]) > 1e-6:
                        print(f"[t={time_now:.4f}] {node.name} mdot_net={mdot_contrib[node]:.6f}")

            # 6. Log States
            for node in processed_nodes:
                node.log_state(time_now)

    # ... (Keep plotting methods the same) ...
    def plot_nodes_overlay(self, nodes, title="Node Comparison", units="SI"):
        fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axs = axs.flatten()
        fig.suptitle(title, fontsize=14)

        for node in nodes:
            time = node.history['time']
            if units == "E":
                axs[0].plot(time, np.array(node.history['P']) / 6894.75729, label=node.name)
                axs[1].plot(time, (np.array(node.history['T']) - 273.15) * 1.8 + 32, label=node.name)
            else:
                axs[0].plot(time, node.history['P'], label=node.name)
                axs[1].plot(time, node.history['T'], label=node.name)
            axs[2].plot(time, node.history['m'], label=node.name)
            axs[3].plot(time, node.history['d'], label=node.name)
            axs[4].plot(time, node.history['Q'], label=node.name)
            axs[5].plot(time, node.history['fill_level'], label=node.name)

        if units == "E":
            axs[0].set_ylabel("Pressure [psi]")
            axs[1].set_ylabel("Temperature [F]")
        else:
            axs[0].set_ylabel("Pressure [Pa]")
            axs[1].set_ylabel("Temperature [K]")
        axs[2].set_ylabel("Mass [kg]")
        axs[3].set_ylabel("Density [kg/m³]")
        axs[4].set_ylabel("Quality [-]")
        axs[4].set_ylim(0, 1)
        axs[5].set_ylabel("Fill level [-]")
        axs[5].set_xlabel("Time [s]")

        for ax in axs:
            ax.legend()
            ax.grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()

    def plot_connections_overlay(self, connections, title="Connection Comparison", units="SI"):
        fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
        axs = axs.flatten()
        fig.suptitle(title, fontsize=14)

        # 1. Unpack Series objects to expose internal components
        items_to_plot = []
        for conn in connections:
            items_to_plot.append((conn, conn.name))
            if hasattr(conn, 'connections'):
                for sub in conn.connections:
                    items_to_plot.append((sub, f"[{conn.name}] {sub.name}"))

        # 2. Plot all items in the flattened list
        for conn, label in items_to_plot:
            time = conn.history['time']
            if not time: continue # Skip if unlogged
            
            if units == "E":
                axs[0].plot(time, conn.history['mdot'], label=label)
                axs[1].plot(time, np.array(conn.history['dP']) / 6894.75729, label=label)
            else:
                axs[0].plot(time, conn.history['mdot'], label=label)
                axs[1].plot(time, conn.history['dP'], label=label)
                
            axs[2].plot(time, np.array(conn.history['CdA']) * 1000000, label=label)
            axs[3].plot(time, conn.history['Hdot'], label=label)
            
            # Filter None values from Quality safely
            Q_clean = [q if q is not None else 0 for q in conn.history['Q']]
            axs[4].plot(time, Q_clean, label=label)
            
            axs[5].plot(time, conn.history['state'], label=label)

        # 3. Format the Graph
        if units == "E":
            axs[0].set_ylabel("mdot [kg/s]")
            axs[1].set_ylabel("dP [psi]")
        else:
            axs[0].set_ylabel("mdot [kg/s]")
            axs[1].set_ylabel("dP [Pa]")
            
        axs[2].set_ylabel("CdA [mm^2]")
        axs[3].set_ylabel("Hdot [J/s]")
        axs[4].set_ylabel("Q [0-1]")
        axs[4].set_ylim(0, 1.1)
        axs[5].set_ylabel("State [-]")
        axs[5].set_xlabel("Time [s]")

        for ax in axs:
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.95])

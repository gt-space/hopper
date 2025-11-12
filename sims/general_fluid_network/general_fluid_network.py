import os
import numpy as np
import matplotlib.pyplot as plt
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import CoolProp.CoolProp as CP

# Initialize REFPROP/CoolProp
try:
    RP = REFPROPFunctionLibrary('C:\\Program Files\\REFPROP') # Modify your install location if necessary
    RP.SETPATHdll(os.environ.get('RPPREFIX', r"C:\Program Files\REFPROP")) # Modify your install location if necessary
    REFPROP = True
    print("Using REFPROP.")
except ValueError:
    REFPROP = False
    print("Using CoolProp.")

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
        return  CP.PropsSI(output, key1, val1, key2, val2, fluid)


class Node():
    """
    Node class. State defined by total density d (kg/m^3) and enthalpy K(J).
    Initialized by fluid, mass m (kg), volume V (L), tempurature T (K), and name.
    """
    def __init__(self, fluid, m, V, T, name="node"):
        self.fluid = fluid
        self.m = float(m) # node mass [kg]
        self.V = float(V) / 1000.0  # convert L -> m^3
        self.name = name

        # Initialize state using T and density computed from m/V
        self.d = self.m / self.V
        # specific enthalpy from (D,T)
        h_spec = PropsSI_auto('H', 'D', self.d, 'T', float(T), self.fluid)  # J/kg
        self.H = self.m * h_spec # total enthalpy [J]
        # derived (will also populate m_l, m_v)
        self._flash_from_DH(self.d, self.H)

        # node history dict initialization
        self.history = {k: [] for k in ["time","Q","P","T","H","h","d","m","m_l","m_v", "fill_level", "s"]}

    def _flash_from_DH(self, d, H):
        """
        Given bulk density d (kg/m3) and total enthalpy H (J),
        compute T, P, h and split m into m_l, m_v if two-phase.
        """
        m = self.m
        if m <= 0:
            # safety floor
            m = 1e-12
            self.m = m

        h = H / m  # specific enthalpy J/kg
        # try to get T and P from (D,H)
        try:
            T = PropsSI_auto('T', 'D', d, 'H', h, self.fluid)
            P = PropsSI_auto('P', 'D', d, 'H', h, self.fluid)
            s = PropsSI_auto('S', 'D', d, 'H', h, self.fluid)
            phase = CP.PhaseSI('D', d, 'H', h, self.fluid) # use only CoolProp here, REFPROP phase lookup behaves weirdly
        except Exception as e:
            raise RuntimeError(f"CoolProp lookup failed in flash: d={d}, h={h}, err={e}") from e

        self.T = T
        self.P = P
        self.h = h
        self.d = d
        self.s = s

        if phase == "twophase":
            Q = PropsSI_auto('Q', 'D', d, 'H', h, self.fluid)  # 0-1
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
            self.m_v = self.m if phase in ("gas", "supercritical") else 0.0
            self.m_l = self.m - self.m_v

            # set phase-specific properties equal to bulk
            self.h_l = self.h_v = self.h
            self.d_l = self.d_v = self.d
            self.fill_level = 1.0 if phase == "liquid" else 0.0

        # safe Cp/Cv/gamma/R in single-phase gas
        try:
            self.Cp = PropsSI_auto('CPMASS', 'D', self.d, 'H', self.h, self.fluid)
            self.Cv = PropsSI_auto('CVMASS', 'D', self.d, 'H', self.h, self.fluid)
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
        self.H += Hdot * dt

        # numerical safety
        if self.m < 1e-12:
            self.m = 1e-12

        # recompute density and flash to get phase split
        d_new = self.m / self.V
        self._flash_from_DH(d_new, self.H)

        # debug print
        # print(f"{self.name}: t-update P={self.P:.1f} Pa, T={self.T:.3f} K, m={self.m:.6f} kg, m_l={self.m_l:.6f}, m_v={self.m_v:.6f}, Q={self.Q}")

    def log_state(self, t=0.0):
        """
        Log node state at each timestep throughout a network sim.
        """
        self.history["time"].append(t)
        self.history["Q"].append(self.Q)
        self.history["P"].append(self.P)
        self.history["T"].append(self.T)
        self.history["H"].append(self.H)
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
        super().__init__(fluid, m=1.0, V=1.0, T=T, name=name)

        # Set fixed ambient conditions
        self.P = P
        self.T = T
        self.h = PropsSI_auto("H", "P", self.P, "T", self.T, fluid)
        self.H = self.h * self.m
        self.d = PropsSI_auto("D", "P", self.P, "T", self.T, fluid)

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


class Accumulator(Node):
    """
    Subclass of Node to represent an accumulator.
    """
    # TODO
    pass


class HeatExchanger(Node):
    """
    Subclass of Node to represnt a heat exchanger.
    """
    # TODO
    pass


class Chamber(Node):
    """
    Subclass of Node to represent a engine combustion chamber.
    """
    # TODO
    pass


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

        # Select phase based on donor fill
        if donor.fill_level > self.location:
            h_stream = donor.h_l
            d_stream = donor.d_l
        else:
            h_stream = donor.h_v
            d_stream = donor.d_v

        abs_dP = abs(dP)
        self.dP = abs_dP  # logging

        donor_phase = CP.PhaseSI('D', donor.d, 'H', donor.h, donor.fluid)

        # --- GAS/CHOKED FLOW ---
        if donor_phase in ("gas", "supercritical") and donor.Cp and donor.Cv and donor.R:
            gamma = donor.gamma
            R = donor.R
            Tdon = donor.T
            crit_factor = ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
            Pcrit = donor.P * crit_factor

            if receiver.P > Pcrit:
                # Unchoked subsonic gas flow
                mdot_mag = self.CdA * donor.P * np.sqrt(2 * abs(1 - (receiver.P / donor.P) ** ((gamma - 1) / gamma)) / (R * Tdon))
            else:
                # Choked
                mdot_mag = self.CdA * donor.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor

            Hdot = mdot_mag * h_stream
        
        # --- LIQUID / TWO-PHASE (Dyer model) ---
        else:
            h_liq = PropsSI_auto('H', 'P', receiver.P, 'Q', 0, donor.fluid)
            h_vap = PropsSI_auto('H', 'P', receiver.P, 'Q', 1, donor.fluid)
            Pv = PropsSI_auto('P', 'T', donor.T, 'Q', 1, donor.fluid)

            # Single-phase incompressible term (SPI)
            mdot_spi = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * abs_dP)

            # Homogeneous equilibrium model term (HEM)
            try:
                h1 = h_stream
                h2 = PropsSI_auto('H', 'P', receiver.P, 'S', donor.s, donor.fluid)
                rho2p = 1.0 / PropsSI_auto('D', 'P', receiver.P, 'Q', 0.5, donor.fluid)
                mdot_hem = self.CdA * rho2p * np.sqrt(2.0 * max(h1 - h2, 1e-9))
            except Exception:
                mdot_hem = mdot_spi

            # Dyer blending factor
            r = 1  # tunable, change based on test data
            kappa = r * (donor.P - receiver.P) / max(Pv - receiver.P, 1e-6) # can also manually set kappa (2 is a good conservative value)

            # Dyer blended mass flow
            mdot_mag = (kappa / (1 + kappa)) * mdot_spi + (1 / (1 + kappa)) * mdot_hem

            # Enthalpy flow rate (assume isenthalpic across connection to find quality)
            # I am aware that an isenthalpic assumption here conflicts with the isentropic assumption for HEM...
            # This is why we blend the models, but generally isenthalpic will be more accurate and conservative...
            # Physically it makes a lot more sense than assuming isentropic (since flashing changes entropy)
            self.Q = PropsSI_auto('Q', 'P', receiver.P, 'H', h_stream, donor.fluid)
            if 0 <= self.Q <= 1:
                Hdot = mdot_mag * (self.Q * h_vap + (1 - self.Q) * h_liq)
            else:
                Hdot = mdot_mag * h_stream

        # Sign convention
        if donor is node1:
            mdot = mdot_mag
        else:
            mdot = -mdot_mag

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
        if upstream.fill_level > self.location:
            h_stream = upstream.h_l
            d_stream = upstream.d_l
        else:
            h_stream = upstream.h_v
            d_stream = upstream.d_v

        donor_phase = CP.PhaseSI('D', upstream.d, 'H', upstream.h, upstream.fluid)
        if donor_phase in ("gas", "supercritical") and upstream.Cp and upstream.Cv and upstream.R:
            gamma = upstream.gamma
            R = upstream.R
            Tdon = upstream.T
            crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
            mdot_mag = self.CdA * upstream.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor
        else:
            mdot_mag = self.CdA * np.sqrt(2.0 * max(d_stream, 1e-6) * effective_dP)

        # Sign convention: positive mdot -> node1 -> node2
        mdot = mdot_mag if upstream is node1 else -mdot_mag

        Hdot = mdot * h_stream + self.qdot
        return mdot, Hdot


class Valve(Connection):
    """
    Subclass of Connection to represent a valve.
    """
    # TODO
    pass


class SharpEdgedOrifice(Connection):
    """
    Subclass of Node to represent a sharp-edged orifice.
    """
    # TODO
    pass


class Engine(Connection):
    """
    Subclass of Connection to represent an engine.
    """    
    # TODO
    pass


class Injector(Connection):
    """
    Subclass of Connection to represent an injector.
    """
    # TODO
    pass


class ThrottleValve(Connection):
    """
    Subclass of Connection to represent a throttle valve.
    """
    def __init__(self, CdA_max, qdot=0, location=0, normal_state=0, checking=True, target_mdot=0, step=0.02, name="throttle_valve"):
        super().__init__(CdA_max*target_mdot*normal_state, qdot, location, normal_state, checking, name)
        self.name = name
        self.target_mdot = target_mdot # target mdot for throttle valve [kg/s]
        self.CdA_max = CdA_max

    def mdot_Hdot(self, node1, node2):
        # TODO ADD DYER
        """
        Return mdot (kg/s), Hdot (J/s) where positive means mass/enthalpy flows node1 -> node2.
        Stream enthalpy is taken from donor's phase at the connection 'location'.
        """
        # check if connection is open
        if not self.state:
            return 0.0, 0.0
    
        dP = node1.P - node2.P
        if self.checking & (dP < 0):
            return 0.0, 0.0
        
        if abs(dP) < 1e-12:
            return 0.0, 0.0

        self.target_mdot = self.state
        # determine donor and receiver
        if dP > 0:
            donor, receiver = node1, node2
        else:
            donor, receiver = node2, node1

        # choose stream phase based on donor.fill_level and connection location
        # if fill_level > location -> stream draws from liquid; otherwise vapor
        if donor.fill_level > self.location:
            # stream is liquid from donor
            h_stream = donor.h_l
            # at receiver pressure
            d_stream = donor.d_l
        else:
            # stream is vapor from donor
            h_stream = donor.h_v
            d_stream = donor.d_v

        abs_dP = abs(dP)
        self.dP = abs_dP # logging

        # allow choked formula only if donor is gas-like (no two-phase)
        donor_phase = CP.PhaseSI('D', donor.d, 'H', donor.h, donor.fluid)
        if donor_phase in ("gas", "supercritical") and donor.Cp and donor.Cv and donor.R:
            # compressible/choked flow estimate
            gamma = donor.gamma
            R = donor.R
            Tdon = donor.T
            crit_factor = ((gamma + 1.0) / 2.0) ** ( - (gamma + 1.0) / (2.0 * (gamma - 1.0)) )
            self.CdA = min(self.CdA_max, self.state / (donor.P / np.sqrt(max(Tdon, 1e-8)) * np.sqrt(gamma / max(R, 1e-12)) * crit_factor))
        else:
            # incompressible-like orifice
            self.CdA = min(self.CdA_max, self.state / np.sqrt(2.0 * max(d_stream, 1e-6) * abs_dP))
            self.Q = PropsSI_auto('Q', 'P', receiver.P, 'H', h_stream, donor.fluid)

            # isentropic assumption test
            # s = PropsSI_auto('S', 'P', donor.P, 'H', h_stream, donor.fluid)
            # self.Q = PropsSI_auto('Q', 'P', receiver.P, 'S', s, donor.fluid)

            # If inside 0–1, it's two-phase
            if 0 <= self.Q <= 1:
                h_liq = PropsSI_auto('H', 'P', receiver.P, 'Q', 0, donor.fluid)
                h_vap = PropsSI_auto('H', 'P', receiver.P, 'Q', 1, donor.fluid)
                m_vap = self.target_mdot * self.Q
                m_liq = self.target_mdot * (1 - self.Q)
                # Hdot adjusted for flashing
                Hdot = m_vap * h_vap + m_liq * h_liq
            else:
                # single-phase at receiver P
                Hdot = self.target_mdot * h_stream

        # sign convention: positive mdot -> node1 -> node2
        if donor is node1:
            mdot = self.target_mdot
        else:
            mdot = - self.target_mdot

        Hdot = mdot * h_stream + self.qdot
        self.mdot = mdot # logging
        self.Hdot = Hdot # logging
        return mdot, Hdot
    

class CheckValve(Connection):
    """
    Subclass of Connection to represent a check valve.
    """
    # TODO
    pass


class Pump(Connection):
    """
    Subclass of Connection to represent a pump.
    """
    # TODO
    pass


class Network():
    """
    Network class. Defined by a graph of connections and nodes representing a fluid network.
    Syntax is as follows: {connection1: (node1, node2), ...}
    """
    def __init__(self, graph):
        self.graph = graph  # {connection: (node1, node2)}

    def sim(self, t, dt, actions={}, verbose_steps=5):
        """
        Runs a transient sim over a fluid network. Takes in total sim runtime, timestep, and an action profile of
        valve set-states and timings.
        Args:
            t: total sim runtime [s] (float)
            dt: timestep duration [s] (float)
            actions: {time [s]: (connection, set_state), ...}, action time must be of the same order of precision as dt.
        """
        steps = int(t / dt)
        for i in range(steps):
            time_now = i * dt
            if time_now in actions:
                actions[time_now][0].state = actions[time_now][1]

            # compute all flows first
            mdot_contrib = {node: 0.0 for _, pair in self.graph.items() for node in pair}
            Hdot_contrib = {node: 0.0 for _, pair in self.graph.items() for node in pair}

            for conn, (n1, n2) in self.graph.items():
                mdot, Hdot = conn.mdot_Hdot(n1, n2)
                # mdot positive => n1 -> n2
                mdot_contrib[n1] -= mdot
                mdot_contrib[n2] += mdot
                Hdot_contrib[n1] -= Hdot
                Hdot_contrib[n2] += Hdot
                conn.log_state(time_now)

            # update nodes simultaneously
            for node in list(mdot_contrib.keys()):
                node.update(mdot_contrib[node], Hdot_contrib[node], dt)
                if i < verbose_steps:
                # print only first few steps to avoid log overload
                    print(f"[t={time_now:.4f}] {node.name} mdot_net={mdot_contrib[node]:.6f}, Hdot_net={Hdot_contrib[node]:.3f}")
                node.log_state(time_now)

    def plot_nodes_overlay(self, nodes, title="Node Comparison", units="SI"):
        """
        Overlay plots of pressure, temperature, mass, density, quality,
        and fill level vs time for all nodes on the same set of subplots.
        Args:
            nodes (list): list of Node objects with histories recorded
            title (str): plot title
            units (str): SI or E
        """
        fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axs = axs.flatten()
        fig.suptitle(title, fontsize=14)

        # Loop over nodes and add each to the plots
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

        # Labels
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

        # Add legends
        for ax in axs:
            ax.legend()
            ax.grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()

    def plot_connections_overlay(self, connections, title="Connection Comparison", units="SI"):
        """
        Overlay plots of CdA, qdot, state, mdot, Hdot
        and dP vs time for all nodes on the same set of subplots.
        Args:
            nodes (list): list of Node objects with histories recorded
            title (str): plot title
            units (str): SI or E
        """
        fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
        axs = axs.flatten()
        fig.suptitle(title, fontsize=14)

        # Loop over nodes and add each to the plots
        for conn in connections:
            time = conn.history['time']
            if units == "E":
                axs[0].plot(time, conn.history['mdot'], label=conn.name)
                axs[1].plot(time, np.array(conn.history['dP']) / 6894.75729, label=conn.name)
            else:
                axs[0].plot(time, conn.history['mdot'], label=conn.name)
                axs[1].plot(time, conn.history['dP'], label=conn.name)
            axs[2].plot(time, conn.history['CdA'], label=conn.name)
            axs[3].plot(time, conn.history['Hdot'], label=conn.name)
            axs[4].plot(time, conn.history['Q'], label=conn.name)
            axs[5].plot(time, conn.history['state'], label=conn.name)

        # Labels
        if units == "E":
            axs[0].set_ylabel("mdot [kg/s]")
            axs[1].set_ylabel("dP [psi]")
        else:
            axs[0].set_ylabel("mdot [kg/s]")
            axs[1].set_ylabel("dP [Pa]")
        axs[2].set_ylabel("CdA [m^2]")
        axs[3].set_ylabel("Hdot [J/s]")
        axs[4].set_ylabel("Q [0-1]")
        axs[5].set_ylabel("State [-]")

        # Add legends
        for ax in axs:
            ax.legend()
            ax.grid(True)

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()
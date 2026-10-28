"""Building blocks for the Hopper P&ID: each element builds native FullFlow parts.

HOPPER_DATA.py describes every flow path as an ordered list of these elements,
from the upstream boundary to the downstream one, e.g.

    LOX_FEED = [
        Line("LOX Tank to OTV Line", inner_diameter=..., length=...),
        ThrottleValve("OTV", max_cda=...),
        Line("OTV to LOX Injector Line", inner_diameter=..., length=...),
        Injector("OINJ", cda=..., length=...),
    ]

Hopper_Prop_System.py walks each list, puts a zero-volume FullFlow junction
(Volume) between consecutive elements, and calls each element's build(). To
add hardware, add an element to the list; nothing else needs to change. The
plots read the per-element tracks written here, so a new element also appears
in the pressure ladders automatically.

Liquid elements (feed lines: tank liquid -> ... -> main chamber)
    Line               FullFlow DarcyWeisbach (with flow inertia) + Churchill friction
    Orifice            FullFlow DischargeCoefficient, fixed CdA (optional passage length = inertia)
    Injector           an Orifice that discharges into the chamber (gives the injector dP / stiffness)
    Valve              FullFlow DischargeCoefficient whose opening follows a time schedule
    ThrottleValve      FullFlow DischargeCoefficient whose CdA a FullFlow Balance sets to hit
                       the thrust-profile flow target (one per feed line)
    CavitatingVenturi  FullFlow CavitatingVenturi (vapour pressure from the propellant)
    Custom             any FullFlow component (or your own Component subclass) you build

Gas elements (pressurant legs: COPV -> ... -> tank ullage)
    GasLine            FullFlow FlowTube (inertia, friction, gas density) + Churchill friction
    Orifice / Valve    fixed or scheduled restriction
    BangBangValve      restriction opened and closed by the tank's bang-bang controller
    CheckValve         restriction that blocks reverse flow into the COPV
    Consecutive gas restrictions carry one flow with no volume between them, so
    they are combined with FullFlow's SeriesCdA into one CompressibleOrifice.

Side paths (vents, reliefs, purges, drains: a named node -> ... -> a sink)
    Branch             restrictions in series from a source node to a sink, one
                       FullFlow CompressibleOrifice (gas) or DischargeCoefficient
                       (liquid) on their SeriesCdA; flow only leaves the source
    ReliefValve        restriction that opens linearly from crack to full-open dP
    Regulator          ideal pressure regulator off a gas node; its outlet is a
                       node that branches can draw from
    Valve(..., schedule=CLOSED) keeps a valve shut for the whole run

Writing your own element: subclass LiquidElement (or use Custom) and implement
build(); it must set self.mass_flow and may set self.cda / self.state for the
ladder plots. FullFlow's own custom-component pattern (a Component subclass with
evaluate_states and balances or dynamics) can be used inside build().
"""

import math

from fullflow import (Balance, CavitatingVenturi as FullFlowCavitatingVenturi, Churchill, DarcyWeisbach,
                      DischargeCoefficient, FlowTube, Sequence, State)


CLOSED = ([0.0], [0.0])   # Valve schedule: shut for the whole run


def tube_id(outer_diameter, wall):
    """Inner diameter [m] of a tube from its outer diameter and wall thickness [m]."""
    return outer_diameter - 2 * wall


def tube_area(inner_diameter):
    """Flow area [m^2] of a round passage."""
    return 3.141592653589793 / 4 * inner_diameter**2


class Tank:
    """A propellant tank: isothermal liquid under a pressurant ullage, with the
    bang-bang controller that holds its pressure.

    propellant          ThermoProp Propellant name (the same name drives the
                        feed-line liquid properties and the combustion map)
    volume [m^3], liquid_mass [kg], liquid_temperature [K]
    pressure [Pa]       starting ullage pressure
    ullage_temperature [K]
    collapse_factor [-] pressurant drawn from the COPV / pressurant that stays
                        in the ullage (ullage collapse: cooling and condensation
                        on the liquid). 1 = no collapse. The COPV loses the full
                        flow; the ullage gains flow / collapse_factor.
    set_pressure, upper_limit, lower_limit [Pa]: the solenoid closes above
                        set + upper and opens below set - lower
    bang_bang_initial_state  1 = solenoid starts open, 0 = closed
    """

    def __init__(self, name, propellant, volume, liquid_mass, liquid_temperature, pressure,
                 ullage_temperature, collapse_factor, set_pressure, upper_limit, lower_limit,
                 bang_bang_initial_state=1.0):
        self.name = name
        self.propellant = propellant
        self.volume = volume
        self.liquid_mass = liquid_mass
        self.liquid_temperature = liquid_temperature
        self.pressure = pressure
        self.ullage_temperature = ullage_temperature
        self.collapse_factor = collapse_factor
        self.set_pressure = set_pressure
        self.upper_limit = upper_limit
        self.lower_limit = lower_limit
        self.bang_bang_initial_state = bang_bang_initial_state


class FeedContext:
    """What an element needs from the network it is built into.

    network       the FullFlow Network
    fluid         FullFlow Lookup of the fluid in this path (density, viscosity, enthalpy, ...)
    target_flow   State with the thrust-profile flow target (liquid feeds only)
    label         "LOX", "Fuel", ... (used in track names)
    """

    def __init__(self, network, fluid, label, target_flow=None):
        self.network = network
        self.fluid = fluid
        self.label = label
        self.target_flow = target_flow


# =========================================================================== #
# Liquid elements
# =========================================================================== #
class LiquidElement:
    """Base class. outlet_node names the junction downstream of this element
    (default "<name> Outlet"); the last element's outlet is the chamber."""

    kind = "liquid"
    closed_at_start = False   # True for elements that start shut (nodes below start at ambient)

    def __init__(self, name, outlet_node=None):
        self.name = name
        self.outlet_node = outlet_node or f"{name} Outlet"
        self.mass_flow = None
        self.upstream_pressure = None
        self.downstream_pressure = None
        self.cda = None       # State/expression [m^2] drawn in the ladder's CdA panel
        self.state = 1.0      # 0-1 shown in the ladder's State panel

    def build(self, context, upstream_pressure, downstream_pressure):
        raise NotImplementedError

    def _connect(self, upstream_pressure, downstream_pressure):
        self.upstream_pressure = upstream_pressure
        self.downstream_pressure = downstream_pressure


class Line(LiquidElement):
    """A liquid tube: FullFlow DarcyWeisbach (dmdot/dt = A/L (dP - friction - gravity))
    with a FullFlow Churchill friction factor, laminar to turbulent.

    inner_diameter [m], length [m], roughness [m], height_change [m] (outlet above inlet).
    """

    def __init__(self, name, inner_diameter, length, roughness=1.5e-5, height_change=0.0, outlet_node=None):
        super().__init__(name, outlet_node)
        self.inner_diameter = inner_diameter
        self.length = length
        self.roughness = roughness
        self.height_change = height_change

    def build(self, context, upstream_pressure, downstream_pressure):
        self._connect(upstream_pressure, downstream_pressure)
        area = tube_area(self.inner_diameter)
        self.mass_flow = State(0.0)
        friction_factor = State(0.02)   # initial guess; Churchill sets it
        Churchill(f"{self.name} Friction", context.network, mass_flow=self.mass_flow,
                  friction_factor=friction_factor, hydraulic_diameter=self.inner_diameter,
                  dynamic_viscosity=context.fluid.dynamic_viscosity, cross_sectional_area=area,
                  roughness=self.roughness)
        DarcyWeisbach(self.name, context.network, mass_flow=self.mass_flow,
                      upstream_pressure=upstream_pressure, downstream_pressure=downstream_pressure,
                      length=self.length, hydraulic_diameter=self.inner_diameter,
                      density=context.fluid.density, cross_sectional_area=area,
                      friction_factor=friction_factor, height_change=self.height_change)
        # Friction-only equivalent CdA, A / sqrt(f L / D) (vehicle_sim.Line's definition).
        self.cda = area / State.maximum(friction_factor * (self.length / self.inner_diameter), 1e-12) ** 0.5


class Orifice(LiquidElement):
    """Fixed restriction: FullFlow DischargeCoefficient, mdot = CdA sqrt(2 rho dP).

    cda [m^2]. length [m] (optional) gives the passage flow inertia, which keeps
    the solver well behaved through zero flow; it does not change the steady dP.
    In a gas leg the same element is a fixed gas restriction.
    """

    def __init__(self, name, cda, length=None, outlet_node=None):
        super().__init__(name, outlet_node)
        self.cda_value = cda
        self.length = length

    def opening(self, context):
        """0-1 multiplier on the CdA (always open for a fixed orifice)."""
        return 1.0

    def build(self, context, upstream_pressure, downstream_pressure):
        self._connect(upstream_pressure, downstream_pressure)
        self.state = self.opening(context)
        component = DischargeCoefficient(
            self.name, context.network, upstream_pressure=upstream_pressure,
            downstream_pressure=downstream_pressure, density=context.fluid.density,
            discharge_coefficient=self.state, cross_sectional_area=self.cda_value,
            length=self.length, mass_flow=State(0.0) if self.length else None)
        self.mass_flow = component.mass_flow
        self.cda = self.cda_value * self.state


class Injector(Orifice):
    """Injector elements: an Orifice whose downstream side is the main chamber.
    Each feed line ends with one; its dP gives the injector stiffness."""


class Valve(Orifice):
    """Valve with a commanded opening: FullFlow DischargeCoefficient whose Cd is
    driven by a FullFlow Sequence.

    cda [m^2] fully open. schedule = (times [s], openings 0-1) on the simulation
    clock, interpolated linearly; None keeps the valve fully open.
    """

    def __init__(self, name, cda, schedule=None, length=None, outlet_node=None):
        super().__init__(name, cda, length, outlet_node)
        self.schedule = schedule
        self.closed_at_start = bool(schedule) and schedule[1][0] == 0.0

    def opening(self, context):
        if not self.schedule:
            return 1.0
        if len(set(self.schedule[1])) == 1:   # constant opening, no Sequence needed
            return float(self.schedule[1][0])
        command = State(float(self.schedule[1][0]))
        Sequence(f"{self.name} Command", context.network, target=command,
                 times=list(self.schedule[0]), values=list(self.schedule[1]))
        return command


class ThrottleValve(LiquidElement):
    """OTV / FTV: a flow controller. FullFlow DischargeCoefficient whose CdA is
    the variable of a FullFlow Balance, so the valve passes the thrust-profile
    flow target (vehicle_sim.py's ThrottleValve "mdot" mode). It cannot open
    past max_cda: a target the feed system cannot deliver shows up as a
    shortfall. Use one per liquid feed.

    max_cda [m^2] fully open. full_scale_mass_flow [kg/s] (optional) only
    scales the ladder's State panel the way the old sim did (target / full scale).
    """

    closed_at_start = True

    def __init__(self, name, max_cda, full_scale_mass_flow=None, outlet_node=None):
        super().__init__(name, outlet_node)
        self.max_cda = max_cda
        self.full_scale_mass_flow = full_scale_mass_flow

    def build(self, context, upstream_pressure, downstream_pressure):
        self._connect(upstream_pressure, downstream_pressure)
        if context.target_flow is None:
            raise ValueError(f"{self.name}: a ThrottleValve needs a flow target (liquid feed lines only).")
        self.area = State(0.0)   # [m^2] solved
        component = DischargeCoefficient(
            self.name, context.network, upstream_pressure=upstream_pressure,
            downstream_pressure=downstream_pressure, density=context.fluid.density,
            discharge_coefficient=1, cross_sectional_area=self.area)
        self.mass_flow = component.mass_flow
        # Residual max((mdot - target) / 0.01 kg/s, (CdA - CdA_max) / 1 mm^2): zero
        # when on target with CdA <= max, or fully open and still short. CdA is
        # unbounded below (targets are never negative; a bound at 0 stalls closing).
        Balance(f"{self.name} Flow Control", context.network, variable=self.area,
                function=State.maximum((self.mass_flow - context.target_flow) / 1.0e-2,
                                       (self.area - self.max_cda) / 1.0e-6))
        self.cda = self.area
        self.target_flow = context.target_flow
        if self.full_scale_mass_flow:
            self.state = State.minimum(State.maximum(context.target_flow / self.full_scale_mass_flow, 0.0), 1.0)


class CavitatingVenturi(LiquidElement):
    """FullFlow CavitatingVenturi: flow set by upstream pressure once the throat
    reaches the propellant's vapour pressure. throat_area [m^2]."""

    def __init__(self, name, throat_area, pressure_recovery_factor=0.85,
                 cavitating_discharge_coefficient=0.94, noncavitating_discharge_coefficient=0.6,
                 outlet_node=None):
        super().__init__(name, outlet_node)
        self.throat_area = throat_area
        self.pressure_recovery_factor = pressure_recovery_factor
        self.cavitating_discharge_coefficient = cavitating_discharge_coefficient
        self.noncavitating_discharge_coefficient = noncavitating_discharge_coefficient

    def build(self, context, upstream_pressure, downstream_pressure):
        self._connect(upstream_pressure, downstream_pressure)
        component = FullFlowCavitatingVenturi(
            self.name, context.network, upstream_pressure=upstream_pressure,
            downstream_pressure=downstream_pressure, density=context.fluid.density,
            throat_area=self.throat_area, vapor_pressure=context.fluid.vapor_pressure,
            pressure_recovery_factor=self.pressure_recovery_factor,
            cavitating_discharge_coefficient=self.cavitating_discharge_coefficient,
            noncavitating_discharge_coefficient=self.noncavitating_discharge_coefficient)
        self.mass_flow = component.mass_flow
        self.cda = self.cavitating_discharge_coefficient * self.throat_area


class Custom(LiquidElement):
    """Any FullFlow component you want in a liquid line.

    build_function(network, name, upstream_pressure, downstream_pressure, fluid)
    must create the component(s) and return the mass-flow State through them.
    fluid is the path's FullFlow Lookup (fluid.density, fluid.dynamic_viscosity, ...).
    Example, a FlowTube with no friction:

        def my_tube(network, name, p_in, p_out, fluid):
            return FlowTube(name, network, mass_flow=State(0.0), upstream_static_pressure=p_in,
                            downstream_static_pressure=p_out, length=0.5, hydraulic_diameter=0.008,
                            upstream_density=fluid.density).mass_flow

        Custom("My Tube", my_tube)
    """

    def __init__(self, name, build_function, outlet_node=None, closed_at_start=False):
        super().__init__(name, outlet_node)
        self.build_function = build_function
        self.closed_at_start = closed_at_start

    def build(self, context, upstream_pressure, downstream_pressure):
        self._connect(upstream_pressure, downstream_pressure)
        self.mass_flow = self.build_function(context.network, self.name, upstream_pressure,
                                             downstream_pressure, context.fluid)
        dp = State.maximum(upstream_pressure - downstream_pressure, 1e-9)
        self.cda = abs(self.mass_flow) / (2.0 * context.fluid.density * dp) ** 0.5


# =========================================================================== #
# Gas elements (pressurant legs)
# =========================================================================== #
class GasLine(LiquidElement):
    """A pressurant tube. inner_diameter [m], length [m], roughness [m].

    inertia=False (default): the tube joins its neighbouring restrictions in the
    leg's FullFlow SeriesCdA as an equivalent CdA, A / sqrt(1 + f L / D), with
    the fully rough friction factor f = (-2 log10(e / 3.7 D))^-2 (vehicle_sim.Line's
    gas-line form). A gas line's acoustic dynamics last milliseconds, far below
    the time step, so this is what the solver can resolve anyway.
    inertia=True: a FullFlow FlowTube (flow inertia, friction from FullFlow's
    Churchill correlation, gas density at both ends) whose gas is stored at its
    outlet node. It resolves pressure waves, but only pays off with a time step
    of a few milliseconds and makes each step much slower.
    """

    kind = "gas_line"

    def __init__(self, name, inner_diameter, length, roughness=1.5e-5, inertia=False, outlet_node=None):
        super().__init__(name, outlet_node)
        self.inner_diameter = inner_diameter
        self.length = length
        self.roughness = roughness
        self.inertia = inertia

    @property
    def cda_value(self):
        """Equivalent CdA [m^2] (fully rough friction), used when inertia=False."""
        friction = (-2.0 * math.log10(self.roughness / (3.7 * self.inner_diameter))) ** -2
        return tube_area(self.inner_diameter) / math.sqrt(1.0 + friction * self.length / self.inner_diameter)

    def opening(self, context):
        return 1.0

    def internal_volume(self):
        """Gas held in the tube [m^3]; stored at its outlet junction."""
        return tube_area(self.inner_diameter) * self.length

    def build_gas(self, network, upstream_pressure, downstream_pressure, upstream_gas, downstream_gas):
        self._connect(upstream_pressure, downstream_pressure)
        area = tube_area(self.inner_diameter)
        self.mass_flow = State(0.0)
        friction_factor = State(0.02)
        Churchill(f"{self.name} Friction", network, mass_flow=self.mass_flow, friction_factor=friction_factor,
                  hydraulic_diameter=self.inner_diameter, dynamic_viscosity=upstream_gas.dynamic_viscosity,
                  cross_sectional_area=area, roughness=self.roughness)
        FlowTube(self.name, network, mass_flow=self.mass_flow, upstream_static_pressure=upstream_pressure,
                 downstream_static_pressure=downstream_pressure, length=self.length,
                 hydraulic_diameter=self.inner_diameter, cross_sectional_area=area,
                 upstream_density=upstream_gas.density, downstream_density=downstream_gas.density,
                 friction_factor=friction_factor)
        self.cda = area / State.maximum(friction_factor * (self.length / self.inner_diameter), 1e-12) ** 0.5


def is_dynamic_line(element):
    """True for a GasLine built as its own FlowTube (inertia=True); every other
    gas element is a restriction lumped into its run's SeriesCdA."""
    return isinstance(element, GasLine) and element.inertia


class BangBangValve(Orifice):
    """Pressurization solenoid opened and closed by its tank's bang-bang
    controller (set pressure and band are given with the tank). cda [m^2]."""

    solenoid = True


class CheckValve(Orifice):
    """Check valve: a restriction that only lets gas flow toward the tank. cda [m^2]."""

    check = True


class ReliefValve(Orifice):
    """Spring relief valve: a restriction whose opening rises linearly from 0 at
    crack_pressure to 1 at full_open_pressure (dP across the valve, i.e. gauge
    pressure when it vents to ambient). cda [m^2] fully open. full_open_pressure
    defaults to 110 % of crack."""

    def __init__(self, name, cda, crack_pressure, full_open_pressure=None, outlet_node=None):
        super().__init__(name, cda, outlet_node=outlet_node)
        self.crack_pressure = crack_pressure
        self.full_open_pressure = full_open_pressure or 1.1 * crack_pressure

    def relief_opening(self, upstream_pressure, downstream_pressure):
        span = self.full_open_pressure - self.crack_pressure
        return ((upstream_pressure - downstream_pressure - self.crack_pressure) / span).clip(0.0, 1.0)


class Regulator:
    """Ideal pressure regulator drawing from a gas node (no droop, no lock-up
    band): its outlet pressure is min(source pressure, ambient + set_pressure)
    at the source gas temperature. Its outlet is the node "<name> Outlet",
    which Branches name as their source; the flow they draw comes out of the
    regulator's source. set_pressure [Pa gauge]. cda [m^2] (optional) is the
    regulator's own flow capacity, put in series with each branch it feeds."""

    def __init__(self, name, source, set_pressure, cda=None):
        self.name = name
        self.source = source
        self.set_pressure = set_pressure
        self.cda_value = cda
        self.outlet_node = f"{name} Outlet"


class Branch:
    """A side path from a named node to a sink: vents, reliefs, purges, drains,
    the torch igniter feeds.

    source  "COPV", "LOX Ullage", "Fuel Ullage" (gas), "LOX Tank Liquid",
            "Fuel Tank Liquid" or any feed junction ("OTV Inlet") (liquid), or
            a Regulator outlet ("PR-REG Outlet")
    sink    "ambient", "Main Chamber", or any feed junction ("OTV Outlet");
            the sink sets the back pressure and the flow leaves the model there
    into_sink  True (liquid branches into "Main Chamber" only): the flow is
            delivered to the chamber and counts in its ox / fuel flow and MR
            (the torch). False: it leaves the model (vents, purges, drains).
    elements  Orifice / Valve / ReliefValve / CheckValve / GasLine in series

    Flow only leaves the source (a branch never back-flows into a tank or the
    COPV). Gas branches are adiabatic restrictions at the source temperature.
    """

    def __init__(self, name, source, elements, sink="ambient", into_sink=False):
        self.name = name
        self.source = source
        self.elements = list(elements)
        self.sink = sink
        self.into_sink = into_sink

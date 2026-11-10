"""Transient Hopper vehicle fluids model (new P&ID), built with FullFlow.

How the network is built
------------------------
FullFlow networks are assembled from components that share State objects: a
component's outlet pressure State is the next component's inlet pressure
State. The flow paths come from HOPPER_DATA.py as ordered element lists
(Hopper_Elements.py); this file walks each list, builds every element's
native FullFlow components, and puts a zero-volume FullFlow Volume (a junction
whose pressure is solved so mass in = mass out) between consecutive elements:

    COPV -> LOX_PRESSURANT_LEG  -> LOX tank ullage
    COPV -> FUEL_PRESSURANT_LEG -> fuel tank ullage
    LOX tank liquid  -> LOX_FEED  -> main chamber -> nozzle
    fuel tank liquid -> FUEL_FEED -> main chamber

Kinds of equations
------------------
  * Storage volumes (COPV, ullages, tank liquids): mass (and energy, for the
    gases) integrated in time with FullFlow's implicit solver.
  * Lines (DarcyWeisbach, FlowTube) and restrictions given a passage length:
    flow inertia, mass flow integrated from the pressure difference.
  * Valves, orifices, nozzle: algebraic flow-vs-pressure relations.
  * Junctions and the main chamber: zero-volume nodes.

Side paths (vents, reliefs, the PR-REG regulator, GN2 purges, OFD and the
torch feeds OIMV / HX-1 / FIMV) are built from BRANCHES / REGULATORS by
build_branches. The torch branches deliver their propellant into the chamber.

Not modelled
------------
Sensors, the torch spark and torch combustion, the ground-side quick
disconnects (PRQD, OQD, FQD), the oxygen diffuser and line priming. The feed
lines start full of liquid, as in vehicle_sim.py. Propellant reaching the
chamber is assumed to ignite at once. Key modelling assumptions are listed in
ASSUMPTIONS.txt.
"""

from functools import lru_cache
import time

from thermoprop import *
from fullflow import *
from fullflow import Balance, CompressibleOrifice, DischargeCoefficient, Lookup, Map, SeriesCdA, State, Volume
from fullflow.System.Components.RocketPerformance import SpecificImpulse   # not re-exported by fullflow
from scipy.constants import *

from HOPPER_DATA import *
from Backend.Hopper_Components import EnginePerformance
from Backend.Hopper_Elements import (BangBangValve, CheckValve, FeedContext, GasLine, Injector, Regen, ReliefValve,
                                     ThrottleValve, Valve, is_dynamic_line)
from scipy.constants import psi
from Backend.Hopper_Control import SampledBangBang, FunctionLookup


# How the throttle-valve flow targets move between thrust-profile points.
#   "step":   hold each point's flow until the next point, moving to the new flow
#             over step_duration() just before it, so the new flow is in place
#             at the point's time (vehicle_sim.py sets each valve once per point).
#   "linear": interpolate the flow linearly between points.
PROFILE_INTERPOLATIONS = ("step", "linear")


def burn_timeline(targets):
    """Key times on the simulation clock [s].

    profile_start  main ignition: OTV/FTV flows reach the first profile point;
                   profile time 0 is here (sequence.main_ignition_s)
    profile_end    last profile point = the shutdown request (OTV starts closing
                   ox_close_after_shutdown_s later, FTV fuel_close_after_shutdown_s later)
    shutdown_end   OTV and FTV both fully closed
    """
    profile_start = PROFILE_START_TIME
    profile_end = profile_start + targets[-1]["time"]
    return profile_start, profile_end, profile_end + max(OX_CLOSE_DELAY, FUEL_CLOSE_DELAY) + MAIN_VALVE_CLOSE_TIME


def valve_window(opened, closed, open_stroke, close_stroke=None):
    """Valve schedule (times, openings): shut, ramp open over open_stroke, hold, ramp shut."""
    close_stroke = open_stroke if close_stroke is None else close_stroke
    open_stroke, close_stroke = max(open_stroke, 1e-3), max(close_stroke, 1e-3)
    return ([0.0, opened, opened + open_stroke, closed, closed + close_stroke], [0.0, 0.0, 1.0, 1.0, 0.0])


def purge_windows(shutdown_request):
    """[(branch name, open time, close time)] on the simulation clock: each purge
    opens its open_after_shutdown_s after the shutdown request (end of the profile)."""
    return [(branch, shutdown_request + delay, shutdown_request + delay + duration)
            for branch, delay, duration in PURGE_SEQUENCE]


def schedule_purges(shutdown_request):
    """Give each purged branch's valve its open window (valves of branches not in the
    sequence keep their own schedule). Returns the time the last purge is shut."""
    by_name = {branch.name: branch for branch in BRANCHES}
    last = shutdown_request
    for name, opened, closed in purge_windows(shutdown_request):
        for element in by_name[name].elements:
            if isinstance(element, Valve) and not isinstance(element, CheckValve):
                element.schedule = valve_window(opened, closed, PURGE_VALVE_STROKE)
        last = max(last, closed + PURGE_VALVE_STROKE)
    return last


def torch_windows(main_purge_end):
    """{valve name: (open, close)} for the torch network on the simulation clock."""
    sequence, ignition = TORCH_SEQUENCE, PROFILE_START_TIME
    close = ignition + sequence["close_after_ignition_s"]
    windows = {}
    for label, start in (("oxidizer", sequence["oxidizer_open_s"]), ("fuel", sequence["fuel_open_s"])):
        for name in TORCH_LEGS[label]["stroke"]:
            windows[name] = (start, close)
    for name in TORCH_LEGS["purge"]["stroke"]:
        windows[name] = (ignition + sequence["purge_open_after_ignition_s"],
                         main_purge_end + sequence["purge_close_after_main_purges_s"])
    return windows


def schedule_torch(main_purge_end):
    """Set the torch valve schedules; returns the time the torch purge is shut."""
    windows, last = torch_windows(main_purge_end), 0.0
    for label in ("oxidizer", "fuel", "purge"):
        leg = TORCH_LEGS[label]
        for element in leg["elements"]:
            if element.name in leg["stroke"]:
                opened, closed = windows[element.name]
                open_stroke, close_stroke = leg["stroke"][element.name]
                element.schedule = valve_window(opened, closed, open_stroke, close_stroke)
                element.closed_at_start = True
                last = max(last, closed + max(close_stroke, 1e-3))
    return last


def step_duration(previous, target):
    """Time [s] OTV/FTV take to move from one profile point to the next: the
    thrust change at THROTTLE_RATE_LIMIT, never less than THROTTLE_STEP_TIME
    and never more than 90 % of the gap between the points."""
    duration = THROTTLE_STEP_TIME
    if THROTTLE_RATE_LIMIT:
        duration = max(duration, abs(target["thrust"] - previous["thrust"]) / THROTTLE_RATE_LIMIT)
    duration = min(duration, 0.9 * (target["time"] - previous["time"]))
    if duration <= 0:
        raise ValueError("step_smoothing_s must be positive in step mode.")
    return duration


def target_schedule(targets, key, interpolation):
    """(times, values) for one commanded quantity on the simulation clock.

    Closed (0) until the valve's start (MAIN_FUEL_OPEN_TIME for fuel,
    MAIN_OX_OPEN_TIME for ox and thrust), ramp to the first profile point by main
    ignition (profile_start), follow the profile, then from the shutdown request
    (profile_end) wait FUEL_CLOSE_DELAY / OX_CLOSE_DELAY and ramp to 0 over
    MAIN_VALVE_CLOSE_TIME. key picks the quantity from the targets
    ("ox_mass_flow", "fuel_mass_flow", "thrust").
    """
    if interpolation not in PROFILE_INTERPOLATIONS:
        raise ValueError(f"interpolation must be one of {PROFILE_INTERPOLATIONS}")

    profile_start, profile_end, shutdown_end = burn_timeline(targets)
    # The opening and closing ramps get a breakpoint every RAMP_BREAKPOINT_SPACING.
    # The ramp is the same straight line, but the solver lands on every breakpoint,
    # so ignition and shutdown are taken in short steps (one 0.05 s step from cold
    # chamber to burning is too big a jump for the Newton solve).
    def ramp(t0, t1, v0, v1):
        n = max(1, int(round((t1 - t0) / RAMP_BREAKPOINT_SPACING)))
        return ([t0 + (t1 - t0) * k / n for k in range(1, n + 1)],
                [v0 + (v1 - v0) * k / n for k in range(1, n + 1)])

    fuel = key == "fuel_mass_flow"
    open_time = MAIN_FUEL_OPEN_TIME if fuel else MAIN_OX_OPEN_TIME
    close_start = profile_end + (FUEL_CLOSE_DELAY if fuel else OX_CLOSE_DELAY)
    times, values = ([0.0, open_time], [0.0, 0.0]) if open_time > 0 else ([0.0], [0.0])
    ramp_times, ramp_values = ramp(open_time, profile_start, 0.0, targets[0][key])
    times += ramp_times
    values += ramp_values

    for previous, target in zip(targets[:-1], targets[1:]):
        point_time = profile_start + target["time"]
        if interpolation == "step":
            times += [point_time - step_duration(previous, target), point_time]
            values += [previous[key], target[key]]
        else:
            times.append(point_time)
            values.append(target[key])

    if close_start > profile_end:
        times.append(close_start)
        values.append(targets[-1][key])
    ramp_times, ramp_values = ramp(close_start, close_start + MAIN_VALVE_CLOSE_TIME, targets[-1][key], 0.0)
    times += ramp_times
    values += ramp_values
    if times[-1] < shutdown_end:
        times.append(shutdown_end)
        values.append(0.0)
    return times, values




# ---- Builders ---- #
LIQUID_PROPERTIES = ("density", "dynamic_viscosity", "enthalpy", "vapor_pressure", "quality")


@lru_cache(maxsize=512)
def _liquid_property_values(propellant, pressure, temperature):
    liquid = Propellant(propellant, pressure=pressure, temperature=temperature)
    return tuple(float(getattr(liquid, name)) for name in LIQUID_PROPERTIES)


def liquid_properties(propellant, pressure, temperature):
    """ThermoProp Propellant liquid properties as plain numbers. Memoised, so
    the Propellant is evaluated once per (propellant, pressure, temperature)
    even though the solver evaluates the network thousands of times."""
    return dict(zip(LIQUID_PROPERTIES, _liquid_property_values(propellant, float(pressure), float(temperature))))


def heated_property_function(propellant, pressure, minimum_temperature, maximum_temperature):
    """Precompute T-dependent properties instead of repeated EOS calls in Newton iterations."""
    import numpy as np
    count = REGEN_PROPERTY_TEMPERATURE_POINTS if maximum_temperature > minimum_temperature else 1
    if count < 1 or (maximum_temperature > minimum_temperature and count < 2):
        raise ValueError("REGEN_PROPERTY_TEMPERATURE_POINTS must be >=2 for a nonzero dT range.")
    temperatures = np.linspace(minimum_temperature, maximum_temperature, count)
    print(f"Building {propellant} liquid-property table: {minimum_temperature:g}–{maximum_temperature:g} K, "
          f"{count} points at {pressure / psi:g} psia...", flush=True)
    start = time.perf_counter()
    rows = []
    last_progress = start
    for index, temperature in enumerate(temperatures):
        rows.append(_liquid_property_values(propellant, pressure, temperature))
        if time.perf_counter() - last_progress >= 5:
            print(f"  Liquid properties: {index + 1}/{count} points ({time.perf_counter() - start:.1f} s)", flush=True)
            last_progress = time.perf_counter()
    values = np.array(rows)
    print(f"Liquid-property table ready in {time.perf_counter() - start:.1f} s.", flush=True)

    def evaluate(temperature):
        return {name: float(np.interp(temperature, temperatures, values[:, i]))
                for i, name in enumerate(LIQUID_PROPERTIES)}
    return evaluate



def track_element(network, element, inlet_enthalpy):
    """Tracks every pressure-ladder quantity of one element, so the plots work
    for any element a user adds."""
    name = element.name
    network.track(f"{name} Mass Flow [kg/s]", element.mass_flow)
    network.track(f"{name} Upstream Pressure [Pa]", element.upstream_pressure)
    network.track(f"{name} Downstream Pressure [Pa]", element.downstream_pressure)
    network.track(f"{name} dP [psi]", (element.upstream_pressure - element.downstream_pressure) / psi)
    network.track(f"{name} CdA [m2]", element.cda)
    network.track(f"{name} State [-]", element.state)
    network.track(f"{name} Inlet Enthalpy [J/kg]", inlet_enthalpy)


def build_tank(network, tank, pressurant):
    """A propellant tank as two FullFlow Volumes sharing the tank volume: an
    isothermal liquid and a pressurant ullage above it (the moving-interface
    pattern of FullFlow's bang-bang example). A Balance moves the liquid volume
    so the liquid surface pressure equals the ullage pressure. Hydrostatic head
    is neglected, as in vehicle_sim.py. Ullage collapse is a collapse factor
    on the pressurant flow, applied by the caller (Vespula's pattern); the
    liquid stays at its fixed temperature.

    pressurant: dict with the pressurant "mass_flow" and "total_enthalpy" into
    the ullage. Returns the liquid Lookup, the liquid and ullage Volumes and
    the ullage gas Lookup.
    """
    # Ullage GN2: starts from (P, T) but is solved in (P, h); near N2's critical
    # point (~492 psia, 126 K) (P, T) states are badly conditioned and (P, h) is not.
    ullage_gas = Lookup(f"{tank.name} Ullage", network, Fluid, PRESSURANT, pressure=tank.pressure,
                        temperature=tank.ullage_temperature, priority=("enthalpy", "temperature"))
    # Liquid propellant properties from the ThermoProp Propellant, evaluated at
    # the tank set pressure and liquid temperature: the liquid is treated as
    # incompressible. The Propellant recomputes a property on every read
    # (0.7 ms for density, ~45 ms for enthalpy), so it is read once here and
    # handed to FullFlow as plain numbers through a function Lookup.
    liquid = FunctionLookup(f"{tank.name} Liquid", network, liquid_properties, propellant=tank.propellant,
                    pressure=tank.set_pressure, temperature=tank.liquid_temperature,
                    outputs=LIQUID_PROPERTIES)

    liquid_volume = State(tank.liquid_mass / liquid.density.value)
    if not 0.0 < liquid_volume.value < tank.volume:
        raise ValueError(f"{tank.name}: liquid inventory must leave positive ullage volume.")
    outflow = State(0.0)   # set by the feed line
    liquid_node = Volume(f"{tank.name} Liquid Node", network, pressure=State(tank.pressure), volume=liquid_volume,
                         density=liquid.density, mass_flow_out=outflow)

    ullage = Volume(f"{tank.name} Ullage Node", network, pressure=ullage_gas.pressure,
                    temperature=ullage_gas.temperature, volume=tank.volume - liquid_node.volume,
                    density=ullage_gas.density, internal_energy=ullage_gas.internal_energy,
                    enthalpy=ullage_gas.enthalpy, energy_variable="enthalpy",
                    mass_flow_in=pressurant["mass_flow"], total_enthalpy_in=pressurant["total_enthalpy"])

    # Scaled to psi so the residual tolerance means the same as elsewhere.
    Balance(f"{tank.name} Liquid Surface Pressure Balance", network, variable=liquid_node.volume,
            function=(liquid_node.pressure - ullage.pressure) / psi)
    return {"liquid": liquid, "liquid_node": liquid_node, "outflow": outflow, "ullage": ullage,
            "ullage_gas": ullage_gas}


def bang_bang_state(network, tank, ullage_pressure, horizon=60.0):
    """The tank's bang-bang solenoid command (1 open, 0 closed), as a FullFlow
    function Sequence: close above set + upper, open below set - lower, else
    keep the last command (hysteresis). Function sequences see the previous
    accepted step, like a flight computer sampling its transducer."""
    command = State(tank.bang_bang_initial_state)
    SampledBangBang(f"{tank.name} Bang Bang", network, ullage_pressure, command,
                   tank.set_pressure, tank.upper_limit, tank.lower_limit, BANG_BANG_SAMPLE_PERIOD, horizon)
    return command


def build_pressurant_leg(network, label, elements, source_gas, source_pressure, receiver_pressure, command):
    """One COPV -> ullage leg from its element list.

    GasLine elements become FullFlow FlowTubes. Each run of neighbouring
    restrictions carries one flow with no volume between them, so it becomes
    one FullFlow CompressibleOrifice whose area is the FullFlow SeriesCdA of
    the run (1/CdA^2 = sum 1/CdA_i^2). The bang-bang command (and any Valve
    schedule) multiplies that orifice's discharge coefficient. Gas nodes
    between segments are zero-volume Volumes with a pressurant Lookup at
    (node pressure, upstream enthalpy): the leg is adiabatic.

    Returns the mass flow drawn from the source, the flow and total enthalpy
    delivered to the receiver, and the junction nodes.
    """
    solenoids = [element for element in elements if isinstance(element, BangBangValve)]
    if len(solenoids) != 1:
        raise ValueError(f"{label} pressurant leg needs exactly one BangBangValve, found {len(solenoids)}.")
    has_check = any(isinstance(element, CheckValve) for element in elements)

    # Split into segments: a GasLine, or a run of restrictions.
    segments = []
    for element in elements:
        if is_dynamic_line(element) or not segments or is_dynamic_line(segments[-1]):
            segments.append(element if is_dynamic_line(element) else [element])
        else:
            segments[-1].append(element)

    upstream_pressure, upstream_gas = source_pressure, source_gas
    nodes, flows = [], []
    for k, segment in enumerate(segments):
        last = k == len(segments) - 1
        if last:
            downstream_pressure, downstream_gas = receiver_pressure, None
        else:
            outlet = segment.outlet_node if is_dynamic_line(segment) else segment[-1].outlet_node
            downstream_pressure = State(PRESSURANT_INITIAL_PRESSURE)
            downstream_gas = Lookup(f"{outlet} Gas", network, Fluid, PRESSURANT, pressure=downstream_pressure,
                                    enthalpy=source_gas.enthalpy)
            nodes.append((outlet, downstream_pressure, downstream_gas))

        if is_dynamic_line(segment):
            if downstream_gas is None:
                raise ValueError(f"{label} pressurant leg must end with a restriction, not a GasLine.")
            segment.build_gas(network, upstream_pressure, downstream_pressure, upstream_gas, downstream_gas)
            track_element(network, segment, upstream_gas.enthalpy)
            flows.append((segment.mass_flow, None))
        else:
            group_name = " + ".join(element.name for element in segment)
            context = FeedContext(network, upstream_gas, label)
            openings = []
            for element in segment:
                if isinstance(element, BangBangValve):
                    opening = element.actuate(context, command, initial=0.0)
                elif isinstance(element, CheckValve):
                    check = FunctionLookup(f"{element.name} Check Command", network,
                                   lambda p1, p2: {"opening": float(p1 > p2)},
                                   p1=upstream_pressure, p2=downstream_pressure, outputs=("opening",))
                    opening = element.actuate(context, check.opening, initial=1.0)
                else:
                    opening = element.opening(context)
                openings.append(opening)
            path = SeriesCdA(f"{group_name} Series CdA", network,
                             effective_areas=[element.cda_value * opening for element, opening in zip(segment, openings)])
            orifice = CompressibleOrifice(
                group_name, network, upstream_total_pressure=upstream_pressure,
                upstream_total_temperature=upstream_gas.temperature, downstream_pressure=downstream_pressure,
                discharge_coefficient=1.0, cross_sectional_area=path.effective_area,
                gas_constant=upstream_gas.gas_constant, specific_heat_ratio=upstream_gas.gamma,
                upstream_static_enthalpy=upstream_gas.enthalpy,
                upstream_static_temperature=upstream_gas.temperature)
            # Check valve: no reverse flow (flow can only go COPV -> tank).
            mass_flow = State.maximum(orifice.mass_flow, 0.0) if has_check else orifice.mass_flow
            network.track(f"{group_name} Effective CdA [m2]", path.effective_area)
            for element, opening in zip(segment, openings):
                element._connect(upstream_pressure, downstream_pressure)
                element.mass_flow = mass_flow
                element.cda = element.cda_value * opening
                element.state = opening
                element.group = group_name
                track_element(network, element, upstream_gas.enthalpy)
            flows.append((mass_flow, orifice.total_enthalpy))

        upstream_pressure, upstream_gas = downstream_pressure, downstream_gas

    # Gas nodes are zero-volume junctions, except next to a line with inertia:
    # there the node stores the line's gas (a FullFlow Volume of the line's
    # internal volume), since an inertial line feeding a solenoid that slams
    # shut needs somewhere to put its moving gas.
    for k, (node_name, node_pressure, node_gas) in enumerate(nodes):
        line = segments[k] if is_dynamic_line(segments[k]) else segments[k + 1]
        storage = {"volume": line.internal_volume(), "density": node_gas.density} if line.inertia else {}
        Volume(node_name, network, pressure=node_pressure, mass_flow_in=flows[k][0],
               mass_flow_out=flows[k + 1][0], **storage)
        network.track(f"{node_name} Pressure [Pa]", node_pressure)

    return {"source_flow": flows[0][0], "mass_flow": flows[-1][0], "total_enthalpy": flows[-1][1],
            "nodes": nodes}


def build_liquid_feed(network, label, elements, tank, chamber_pressure, target_flow):
    """One tank -> chamber feed from its element list: each element builds its
    own FullFlow components; a zero-volume Volume joins each pair. Nodes behind
    a closed element (the throttle valve at t = 0) start at ambient, the rest at
    tank pressure. Returns the tank outflow, the chamber inflow and the nodes."""
    throttles = [element for element in elements if isinstance(element, ThrottleValve)]
    if len(throttles) != 1:
        raise ValueError(f"{label} feed needs exactly one ThrottleValve, found {len(throttles)}.")
    if not isinstance(elements[-1], Injector):
        raise ValueError(f"{label} feed must end with an Injector (it discharges into the chamber).")

    context = FeedContext(network, tank["liquid"], label, target_flow)
    liquid_temperature = State(LOX_TANK.liquid_temperature if label == "LOX" else FUEL_TANK.liquid_temperature)
    temperature_min = temperature_max = liquid_temperature.value
    upstream_pressure, closed, nodes = tank["liquid_node"].pressure, False, []
    for i, element in enumerate(elements):
        closed = closed or element.closed_at_start
        if i == len(elements) - 1:
            downstream_pressure = chamber_pressure
        else:
            downstream_pressure = State(AMBIENT_PRESSURE if closed else tank["liquid_node"].pressure.value)
            nodes.append((element.outlet_node, downstream_pressure))
        element.build(context, upstream_pressure, downstream_pressure)
        track_element(network, element, context.fluid.enthalpy)
        network.track(f"{element.name} Inlet Temperature [K]", liquid_temperature)
        if isinstance(element, Regen):
            inlet_temperature, inlet_enthalpy = liquid_temperature, context.fluid.enthalpy
            temperature = FunctionLookup(f"{element.name} Outlet Temperature", network,
                                 lambda inlet, pc, part=element: {"temperature": part.outlet_temperature(inlet, pc)},
                                 inlet=liquid_temperature, pc=chamber_pressure, outputs=("temperature",))
            liquid_temperature = temperature.temperature
            spec = LOX_TANK if label == "LOX" else FUEL_TANK
            temperature_min += min(dt for _, dt in element.delta_t_curve)
            temperature_max += max(dt for _, dt in element.delta_t_curve)
            properties = heated_property_function(spec.propellant, spec.set_pressure, temperature_min, temperature_max)
            context.fluid = FunctionLookup(f"{element.name} Heated Liquid", network, properties,
                                   temperature=liquid_temperature, outputs=LIQUID_PROPERTIES)
            network.track(f"{element.name} dT [K]", liquid_temperature - inlet_temperature)
            network.track(f"{element.name} Heat Pickup [W]", element.mass_flow * (context.fluid.enthalpy - inlet_enthalpy))
        network.track(f"{element.name} Outlet Temperature [K]", liquid_temperature)
        if i < len(elements) - 1:
            network.track(f"{element.outlet_node} Temperature [K]", liquid_temperature)
            network.track(f"{element.outlet_node} Density [kg/m3]", context.fluid.density)
        upstream_pressure = downstream_pressure

    junctions = {}
    for i, (node_name, node_pressure) in enumerate(nodes):
        junction = Volume(node_name, network, pressure=node_pressure, mass_flow_in=elements[i].mass_flow,
                          mass_flow_out=elements[i + 1].mass_flow)
        junctions[node_name] = (junction, elements[i + 1].mass_flow)   # side branches add to the outflow
        network.track(f"{node_name} Pressure [Pa]", node_pressure)

    throttle = throttles[0]
    network.track(f"{throttle.name} CdA [mm2]", throttle.cda / milli**2)
    network.track(f"{throttle.name} Target Mass Flow [kg/s]", target_flow)
    network.track(f"{label} Injector dP [psi]",
                  (elements[-1].upstream_pressure - chamber_pressure) / psi)
    return {"tank_outflow": elements[0].mass_flow, "chamber_inflow": elements[-1].mass_flow, "nodes": nodes,
            "junctions": junctions, "outlet_temperature": liquid_temperature}


def regulator_outlet_volume(regulator, branches, minimum=1.0e-6):
    """Gas volume [m^3] always connected to a regulator outlet: the GasLines of every
    branch it feeds, up to that branch's first valve (the purge lines run full-length
    from the regulator to their solenoids). At least minimum (1 mL)."""
    volume = 0.0
    for branch in branches:
        if branch.source != regulator.outlet_node:
            continue
        for element in branch.elements:
            if isinstance(element, GasLine):
                volume += element.internal_volume()
            else:
                break
    return max(volume, minimum)


def build_regulators(network, regulators, nodes, outflows, branches=()):
    """Regulator outlet nodes. Adds "<name> Outlet" to nodes; the regulator's draw
    is appended to outflows[source store].

    cda None: ideal, outlet = min(source, target), at the source gas state.
    cda set: the regulator is a FullFlow CompressibleOrifice (source -> outlet)
    whose opening rises linearly from 0 at the target outlet pressure to 1 at
    droop below it (a proportional regulator); with the valve fully open the outlet
    falls further (capacity limited).
    The outlet node is shared by every branch the regulator feeds (finish_regulators
    sets its outflow) and stores the gas in the lines between the regulator and their
    valves (regulator_outlet_volume); its gas is the source gas throttled at constant
    enthalpy. (An exact-set regulator on a zero-volume outlet stalls the solver when a
    purge opens or closes on a fully open regulator.) The target is the set plus the supply-pressure
    effect, set + effect * (reference inlet - inlet)."""
    for regulator in regulators:
        source = nodes[regulator.source]
        if not source["gas"]:
            raise ValueError(f"{regulator.name}: a Regulator draws from a gas node, not {regulator.source}.")
        inlet, fluid = source["pressure"], source["fluid"]
        reference = float(inlet.value) if hasattr(inlet, "value") else float(inlet)
        target = regulator.outlet_target(inlet, AMBIENT_PRESSURE, reference)
        name = regulator.outlet_node
        if not regulator.cda_value:
            outlet_pressure = State.minimum(inlet, target)
            nodes[name] = dict(source, pressure=outlet_pressure, regulator=regulator)
            network.track(f"{name} Pressure [psia]", outlet_pressure / psi)
            continue
        start = min(reference, AMBIENT_PRESSURE + regulator.set_pressure)
        outlet_pressure = State(start)
        gas = Lookup(f"{name} Gas", network, Fluid, PRESSURANT, pressure=outlet_pressure, enthalpy=fluid.enthalpy)
        # Proportional (droop) regulator: shut at the target, fully open droop below it.
        opening = ((target - outlet_pressure) / regulator.droop).clip(0.0, 1.0)
        valve = CompressibleOrifice(
            regulator.name, network, upstream_total_pressure=inlet, upstream_total_temperature=fluid.temperature,
            downstream_pressure=outlet_pressure, discharge_coefficient=1.0,
            cross_sectional_area=regulator.cda_value * opening, gas_constant=fluid.gas_constant,
            specific_heat_ratio=fluid.gamma, upstream_static_enthalpy=fluid.enthalpy,
            upstream_static_temperature=fluid.temperature)
        junction = Volume(name, network, pressure=outlet_pressure, mass_flow_in=valve.mass_flow,
                          volume=regulator_outlet_volume(regulator, branches), density=gas.density)
        nodes[name] = {"pressure": outlet_pressure, "fluid": gas, "gas": True, "store": name, "volume": junction,
                       "label": source.get("label")}
        outflows.setdefault(source["store"], []).append(valve.mass_flow)
        network.track(f"{name} Pressure [psia]", outlet_pressure / psi)
        network.track(f"{name} Target Pressure [psia]", target / psi)
        network.track(f"{name} Temperature [K]", gas.temperature)
        network.track(f"{regulator.name} Opening", opening)
        network.track(f"{regulator.name} Mass Flow [kg/s]", valve.mass_flow)


def finish_regulators(regulators, nodes, outflows):
    """Every flow drawn from a capacity-limited regulator's outlet leaves its junction."""
    for regulator in regulators:
        node = nodes[regulator.outlet_node]
        if "volume" in node:
            flows = outflows.get(regulator.outlet_node, [])
            node["volume"].mass_flow_out = sum(flows[1:], flows[0]) if flows else State(0.0)


def build_branches(network, regulators, branches, nodes, sinks, outflows=None, finish=True):
    """Regulators and side branches (vents, reliefs, purges, drains) from
    HOPPER_DATA.py. finish=False leaves the regulator outlets open for more
    draws (the torch purge); call finish_regulators afterwards.

    nodes: {name: {"pressure", "fluid", "gas" (True/False), "store", "label"}}
    for the gas and liquid nodes a branch can draw from; "store" names the
    Volume that loses the flow, "label" ("LOX" / "Fuel") the propellant. sinks: {name: pressure}. Each branch's restrictions are
    one FullFlow SeriesCdA feeding one CompressibleOrifice (gas) or
    DischargeCoefficient (liquid); Valve schedules and ReliefValve openings
    multiply its discharge coefficient, and the flow is clipped at zero so it
    only leaves the source. Returns {store name: [outflow, ...]} and the
    flows delivered into the chamber, {"LOX": [...], "Fuel": [...]}.
    """
    outflows = {} if outflows is None else outflows
    delivered = {}
    build_regulators(network, regulators, nodes, outflows, branches)

    for branch in branches:
        if branch.source not in nodes:
            raise ValueError(f"{branch.name}: unknown source {branch.source!r}; one of {sorted(nodes)}.")
        if branch.sink not in sinks:
            raise ValueError(f"{branch.name}: unknown sink {branch.sink!r}; one of {sorted(sinks)}.")
        if any(is_dynamic_line(element) for element in branch.elements):
            raise ValueError(f"{branch.name}: GasLine(inertia=True) is not supported in a Branch.")
        source = nodes[branch.source]
        upstream, downstream, fluid = source["pressure"], sinks[branch.sink], source["fluid"]

        context = FeedContext(network, fluid, branch.name)
        openings = [element.relief_opening(upstream, downstream) if isinstance(element, ReliefValve)
                    else element.opening(context) for element in branch.elements]
        discharge_coefficient = 1.0
        for opening in openings:
            discharge_coefficient = discharge_coefficient * opening
        areas = [element.cda_value for element in branch.elements]
        regulator = source.get("regulator")
        if regulator is not None and regulator.cda_value:
            areas.insert(0, regulator.cda_value)
        path = SeriesCdA(f"{branch.name} Series CdA", network, effective_areas=areas)

        if source["gas"]:
            component = CompressibleOrifice(
                branch.name, network, upstream_total_pressure=upstream, upstream_total_temperature=fluid.temperature,
                downstream_pressure=downstream, discharge_coefficient=discharge_coefficient,
                cross_sectional_area=path.effective_area, gas_constant=fluid.gas_constant,
                specific_heat_ratio=fluid.gamma)
        else:
            component = DischargeCoefficient(
                branch.name, network, upstream_pressure=upstream, downstream_pressure=downstream,
                density=fluid.density, discharge_coefficient=discharge_coefficient,
                cross_sectional_area=path.effective_area)
        mass_flow = State.maximum(component.mass_flow, 0.0)
        outflows.setdefault(source["store"], []).append(mass_flow)
        if branch.into_sink:
            if source["gas"] or branch.sink != "Main Chamber":
                raise ValueError(f"{branch.name}: into_sink is only for liquid branches into the Main Chamber.")
            delivered.setdefault(source["label"], []).append(mass_flow)

        network.track(f"{branch.name} Mass Flow [kg/s]", mass_flow)
        network.track(f"{branch.name} Effective CdA [m2]", path.effective_area * discharge_coefficient)
        for element, opening in zip(branch.elements, openings):
            element._connect(upstream, downstream)
            element.mass_flow = mass_flow
            element.cda = element.cda_value
            element.state = opening
            element.group = branch.name
            track_element(network, element, fluid.enthalpy)
    if finish:
        finish_regulators(regulators, nodes, outflows)
    return outflows, delivered


def _torch_path(network, name, elements, shared, upstream, fluid, downstream, gas):
    """Restrictions in series (SeriesCdA) from a node into the torch chamber;
    scheduled valves scale the discharge coefficient; no reverse flow. shared:
    elements of another leg the flow also passes (counted in the CdA, not tracked)."""
    context = FeedContext(network, fluid, name)
    openings = [element.opening(context) for element in elements]
    discharge = 1.0
    for opening in openings:
        discharge = discharge * opening
    path = SeriesCdA(f"{name} Series CdA", network,
                     effective_areas=[element.cda_value for element in list(elements) + list(shared)])
    if gas:
        component = CompressibleOrifice(
            name, network, upstream_total_pressure=upstream, upstream_total_temperature=fluid.temperature,
            downstream_pressure=downstream, discharge_coefficient=discharge, cross_sectional_area=path.effective_area,
            gas_constant=fluid.gas_constant, specific_heat_ratio=fluid.gamma)
    else:   # the flow is assigned before the first solve: the chamber MR reads it
        component = DischargeCoefficient(
            name, network, upstream_pressure=upstream, downstream_pressure=downstream, density=fluid.density,
            discharge_coefficient=discharge, cross_sectional_area=path.effective_area, mass_flow=State(0.0))
    mass_flow = State.maximum(component.mass_flow, 0.0)
    network.track(f"{name} Mass Flow [kg/s]", mass_flow)
    network.track(f"{name} Effective CdA [m2]", path.effective_area * discharge)
    for element, opening in zip(elements, openings):
        element._connect(upstream, downstream)
        element.mass_flow, element.cda, element.state, element.group = mass_flow, element.cda_value, opening, name
        track_element(network, element, fluid.enthalpy)
    return mass_flow


def build_torch(network, feeds, tanks, nodes, chamber_pressure, outflows):
    """The torch igniter from TORCH_LEGS: LOX from a LOX feed node and fuel from a fuel
    feed node into a zero-volume torch chamber that burns on the main chamber map;
    the TFPV purge gas (from its regulator outlet, through the shared torch fuel
    injector) also enters it. The torch throat carries the mixture into the main
    chamber (ideal, subsonic or choked, either direction). Draws are appended to
    outflows (feed junctions, regulator outlet)."""
    from Backend.Hopper_Components import torch_throat_flow
    legs = TORCH_LEGS
    torch_pressure = State(AMBIENT_PRESSURE)
    flows = {}
    for label, side, title in (("oxidizer", "LOX", "Torch Ox"), ("fuel", "Fuel", "Torch Fuel")):
        leg = legs[label]
        feed_nodes = dict(feeds[side]["nodes"])
        if leg["source"] not in feed_nodes:
            raise ValueError(f"torch {label} leg: source {leg['source']!r} must be a {side} feed junction.")
        flows[side] = _torch_path(network, title, leg["elements"], [], feed_nodes[leg["source"]],
                                  tanks[side]["liquid"], torch_pressure, gas=False)
        outflows.setdefault(leg["source"], []).append(flows[side])
    purge = legs["purge"]
    source = nodes[purge["source"]]
    purge_flow = _torch_path(network, "Torch Purge", purge["elements"], purge["shared"], source["pressure"],
                             source["fluid"], torch_pressure, gas=True)
    outflows.setdefault(source["store"], []).append(purge_flow)

    floor = MIXTURE_RATIO_FLOW_FLOOR
    torch_mr = ((flows["LOX"] + MIXTURE_RATIO * floor) / (flows["Fuel"] + floor)).clip(
        CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT, CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT)
    torch_gas = Map.from_hdf5("Torch Gas", network, CHAMBER_MAP_FILENAME, group="chamber",
                              inputs={"chamber_pressure": State.maximum(torch_pressure, CHAMBER_MAP_PRESSURE_LOWER_LIMIT),
                                      "mixture_ratio": torch_mr}, extrapolate=True)
    purge_gas = source["fluid"]
    throat = FunctionLookup("Torch Throat", network, torch_throat_flow, p_up=torch_pressure, p_down=chamber_pressure,
                            combustion_flow=flows["LOX"] + flows["Fuel"], purge_flow=purge_flow,
                            gamma_c=torch_gas.gamma, r_c=torch_gas.gas_constant, t_c=torch_gas.temperature,
                            gamma_n=purge_gas.gamma, r_n=purge_gas.gas_constant, t_n=purge_gas.temperature,
                            area=legs["throat_area"] * legs["throat_discharge_coefficient"],
                            band=NOZZLE_LINEAR_PRESSURE_DROP, outputs=("mass_flow", "temperature"))
    Volume(legs["chamber_name"], network, pressure=torch_pressure,
           mass_flow_in=flows["LOX"] + flows["Fuel"] + purge_flow, mass_flow_out=throat.mass_flow)
    network.track("Torch Chamber Pressure [psia]", torch_pressure / psi)
    network.track("Torch Over Main Pressure [psi]", (torch_pressure - chamber_pressure) / psi)
    network.track("Torch Mixture Ratio", torch_mr)
    network.track("Torch Gas Temperature [K]", throat.temperature)
    network.track("Torch Throat Mass Flow [kg/s]", throat.mass_flow)
    network.track("Torch LOX Mass Flow [kg/s]", flows["LOX"])
    network.track("Torch Fuel Mass Flow [kg/s]", flows["Fuel"])
    return {"pressure": torch_pressure, "ox": flows["LOX"], "fuel": flows["Fuel"], "purge": purge_flow,
            "outflow": throat.mass_flow}


def build_lumped_copv(network):
    """All bottles as one GN2 volume: identical bottles at the same state blow down
    together, so one volume gives the same answer without the stiff header. Each
    bottle's tracks are this volume's state and its share of the gas."""
    if len({(b.pressure, b.temperature) for b in COPVS}) != 1:
        raise ValueError('copv_model "lumped" needs every bottle at the same pressure and temperature; '
                         'use "header" for bottles that differ.')
    gas = Lookup("COPV Pressurant Gas", network, Fluid, PRESSURANT,
                 pressure=COPVS[0].pressure, temperature=COPVS[0].temperature)
    copv = Volume("COPV", network, pressure=gas.pressure, temperature=gas.temperature,
                  volume=COPV_VOLUME, density=gas.density, internal_energy=gas.internal_energy,
                  enthalpy=gas.enthalpy, energy_variable="temperature")
    for spec in COPVS:
        network.track(f"{spec.name} Pressure [psia]", copv.pressure / psi)
        network.track(f"{spec.name} Temperature [K]", copv.temperature)
        network.track(f"{spec.name} Gas Mass [kg]", copv.mass * (spec.volume / COPV_VOLUME))
    return copv, gas, copv.mass


def build_copvs(network):
    """Independent bottles feeding an adiabatic common header with gas storage."""
    if not COPVS or len({b.name for b in COPVS}) != len(COPVS):
        raise ValueError("COPVS must contain at least one bottle with unique names.")
    if COPV_MODEL == "lumped":
        return build_lumped_copv(network)
    pressure = State(min(b.pressure for b in COPVS))
    bottles, flows = [], []
    enthalpy_flow = 0.0
    for spec in COPVS:
        gas = Lookup(f"{spec.name} Gas", network, Fluid, PRESSURANT,
                     pressure=spec.pressure, temperature=spec.temperature)
        bottle = Volume(spec.name, network, pressure=gas.pressure, temperature=gas.temperature,
                        volume=spec.volume, density=gas.density, internal_energy=gas.internal_energy,
                        enthalpy=gas.enthalpy, energy_variable="temperature")
        outlet = CompressibleOrifice(f"{spec.name} Header Connection", network,
                    upstream_total_pressure=bottle.pressure, upstream_total_temperature=gas.temperature,
                    downstream_pressure=pressure, cross_sectional_area=spec.outlet_cda,
                    discharge_coefficient=1.0, gas_constant=gas.gas_constant, specific_heat_ratio=gas.gamma,
                    upstream_static_enthalpy=gas.enthalpy, upstream_static_temperature=gas.temperature)
        flow = State.maximum(outlet.mass_flow, 0.0)
        bottle.mass_flow_out = flow
        flows.append(flow)
        enthalpy_flow = enthalpy_flow + flow * gas.enthalpy
        bottles.append((spec, bottle, gas))
        network.track(f"{spec.name} Pressure [psia]", bottle.pressure / psi)
        network.track(f"{spec.name} Temperature [K]", bottle.temperature)
        network.track(f"{spec.name} Gas Mass [kg]", bottle.mass)
        network.track(f"{spec.name} Outflow [kg/s]", flow)
    total_flow = sum(flows)
    if COPV_HEADER_VOLUME <= 0:
        raise ValueError("COPV_HEADER_VOLUME must be positive so trapped header gas has a physical state.")
    gas = Lookup("COPV Pressurant Gas", network, Fluid, PRESSURANT, pressure=pressure,
                 temperature=sum(b.temperature for b in COPVS) / len(COPVS))
    # FullFlow takes total_enthalpy_in per kg (energy in = mass_flow_in * total_enthalpy_in),
    # so the header gets the flow-weighted mix of the bottle enthalpies, not the enthalpy flow.
    header = Volume("COPV Header", network, pressure=pressure, temperature=gas.temperature,
                    volume=COPV_HEADER_VOLUME, density=gas.density, internal_energy=gas.internal_energy,
                    enthalpy=gas.enthalpy, energy_variable="temperature", mass_flow_in=total_flow,
                    total_enthalpy_in=enthalpy_flow / State.maximum(total_flow, 1.0e-9))
    network.track("COPV Header Gas Mass [kg]", header.mass)
    return header, gas, sum(b.mass for _, b, _ in bottles)




def HopperPropSystem(targets: list,
                     dt: float = PROP_DT,
                     t_final: float = None,
                     filename: str = PROP_FILENAME,
                     chamber_map_filename: str = CHAMBER_MAP_FILENAME,
                     profile_interpolation: str = PROP_PROFILE_INTERPOLATION,
                     rtol: float = PROP_RTOL,
                     jacobian_method: str = PROP_JACOBIAN_METHOD,
                     save_dt: float = PROP_SAVE_DT,
                     verbose: bool = False,
                     statistics: bool = False):
    """Build the Hopper prop system network from HOPPER_DATA.py and run the transient.

    targets                output of HopperThrustTargets (flows per profile point)
    dt                     solver time step [s]
    t_final                end time [s]; default is shutdown + TAIL_OFF_TIME
    filename               FullFlow HDF5 output
    profile_interpolation  "step" (as vehicle_sim.py) or "linear", see above
    rtol                   per-step residual acceptance tolerance
    jacobian_method        finite-difference Jacobian ("2-point" is ~2x faster than "3-point")
    save_dt                output interval [s] (None saves every step)
    """
    profile_start, profile_end, shutdown_end = burn_timeline(targets)
    # The purges open at their delays after the shutdown request (purge section of the config).
    purge_end = schedule_purges(profile_end) if PURGE_IN_TRANSIENT and PURGE_SEQUENCE else shutdown_end
    # Torch startup and the TFPV purge, which stays open until after the main purges.
    torch_end = schedule_torch(purge_end) if TORCH_LEGS else 0.0
    if t_final is None:
        t_final = max(shutdown_end, purge_end, torch_end) + TAIL_OFF_TIME

    build_start = time.perf_counter()
    print("Building flow network and fluid properties...", flush=True)
    Hopper = Network("Hopper Prop System")


    # ---- Commands ---- #
    # Throttle-valve flow targets on the simulation clock (see target_schedule).
    # THRUST_TARGET is only tracked for plotting; it does not drive anything.
    LOX_TARGET_FLOW, FUEL_TARGET_FLOW, THRUST_TARGET = State(0.0), State(0.0), State(0.0)
    for target_state, key, name in [(LOX_TARGET_FLOW, "ox_mass_flow", "OTV Flow Target"),
                                    (FUEL_TARGET_FLOW, "fuel_mass_flow", "FTV Flow Target"),
                                    (THRUST_TARGET, "thrust", "Thrust Target")]:
        times, values = target_schedule(targets, key, profile_interpolation)
        Sequence(name, Hopper, target=target_state, times=times, values=values)


    # ---- Pressurant: COPVs ---- #
    # One lumped GN2 volume, or bottles feeding a common header (copv_model); mass
    # and energy integrated, no wall heat transfer, so this is an adiabatic blowdown.
    COPV, COPVGas, COPV_GAS_MASS = build_copvs(Hopper)


    # ---- Tanks and pressurant legs ---- #
    tanks, legs = {}, {}
    for label, tank, leg in [("LOX", LOX_TANK, LOX_PRESSURANT_LEG), ("Fuel", FUEL_TANK, FUEL_PRESSURANT_LEG)]:
        pressurant = {"mass_flow": State(0.0), "total_enthalpy": State(0.0)}   # filled in below
        tanks[label] = build_tank(Hopper, tank, pressurant)
        command = bang_bang_state(Hopper, tank, tanks[label]["ullage_gas"].pressure, t_final)
        legs[label] = build_pressurant_leg(Hopper, label, leg, COPVGas, COPV.pressure,
                                           tanks[label]["ullage"].pressure, command)
        # Ullage collapse: the COPV supplies the full leg flow, the ullage keeps 1/CF of it.
        tanks[label]["ullage"].mass_flow_in = legs[label]["mass_flow"] / tank.collapse_factor
        tanks[label]["ullage"].total_enthalpy_in = legs[label]["total_enthalpy"]
        legs[label]["command"] = command


    # ---- Feed lines ---- #
    PC = State(AMBIENT_PRESSURE)   # [Pa] main chamber pressure (solved)
    feeds = {
        "LOX": build_liquid_feed(Hopper, "LOX", LOX_FEED, tanks["LOX"], PC, LOX_TARGET_FLOW),
        "Fuel": build_liquid_feed(Hopper, "Fuel", FUEL_FEED, tanks["Fuel"], PC, FUEL_TARGET_FLOW),
    }


    # ---- Vents, reliefs, purges, drains ---- #
    nodes = {"COPV": {"pressure": COPV.pressure, "fluid": COPVGas, "gas": True, "store": "COPV"}}
    sinks = {"ambient": AMBIENT_PRESSURE, "Main Chamber": PC}
    for label in tanks:
        nodes[f"{label} Ullage"] = {"pressure": tanks[label]["ullage_gas"].pressure,
                                    "fluid": tanks[label]["ullage_gas"], "gas": True, "store": f"{label} Ullage"}
        nodes[f"{label} Tank Liquid"] = {"pressure": tanks[label]["liquid_node"].pressure,
                                         "fluid": tanks[label]["liquid"], "gas": False, "store": f"{label} Liquid",
                                         "label": label}
        for name, pressure in feeds[label]["nodes"]:
            sinks[name] = pressure
            nodes[name] = {"pressure": pressure, "fluid": tanks[label]["liquid"], "gas": False, "store": name,
                           "label": label}
    branch_flows, delivered = build_branches(Hopper, REGULATORS, BRANCHES, nodes, sinks, finish=False)
    torch = build_torch(Hopper, feeds, tanks, nodes, PC, branch_flows) if TORCH_LEGS else None
    finish_regulators(REGULATORS, nodes, branch_flows)

    def total(flows):
        return sum(flows[1:], flows[0]) if flows else 0.0

    COPV.mass_flow_out = total([legs["LOX"]["source_flow"], legs["Fuel"]["source_flow"]]
                               + branch_flows.get("COPV", []))
    for label in feeds:
        tanks[label]["liquid_node"].mass_flow_out = total([feeds[label]["tank_outflow"]]
                                                          + branch_flows.get(f"{label} Liquid", []))
        if branch_flows.get(f"{label} Ullage"):
            tanks[label]["ullage"].mass_flow_out = total(branch_flows[f"{label} Ullage"])
        for name, (junction, next_flow) in feeds[label]["junctions"].items():
            if branch_flows.get(name):
                junction.mass_flow_out = total([next_flow] + branch_flows[name])


    # ---- Main chamber ---- #
    # Zero-volume node, as in Vespula: Pc is whatever makes the nozzle outflow
    # equal the injector inflow. (vehicle_sim.Engine instead relaxed Pc 30 %
    # per 0.1 s step toward that value.)
    # Main injectors plus anything a branch delivers into the chamber (the torch).
    # The torch's LOX and fuel count in the chamber's flow and MR (they reach it as torch
    # products through the torch throat). The TFPV nitrogen pressurizes the torch and
    # flows out of its throat, but is not counted in the main chamber: the chamber map
    # would treat it as hot combustion gas and overstate Pc and thrust (~1.8 %).
    torch_ox, torch_fuel = ([torch["ox"]], [torch["fuel"]]) if torch else ([], [])
    OX_INTO_CHAMBER = total([feeds["LOX"]["chamber_inflow"]] + delivered.get("LOX", []) + torch_ox)
    FUEL_INTO_CHAMBER = total([feeds["Fuel"]["chamber_inflow"]] + delivered.get("Fuel", []) + torch_fuel)
    # Chamber MR from the ox and fuel entering, kept defined at zero flow (see
    # MIXTURE_RATIO_FLOW_FLOOR) and clipped to the map range.
    MR = ((OX_INTO_CHAMBER + MIXTURE_RATIO * MIXTURE_RATIO_FLOW_FLOOR)
          / (FUEL_INTO_CHAMBER + MIXTURE_RATIO_FLOW_FLOOR)).clip(
        CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT, CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT)

    # Gas map, nozzle, efficiencies and thrust (shared with the targeter).
    Engine = EnginePerformance("Main Chamber", Hopper, chamber_pressure=PC, mixture_ratio=MR, sampled_efficiency=True,
                               chamber_map_filename=chamber_map_filename)
    Chamber = Volume("Main Chamber", Hopper, pressure=PC, mass_flow_in=OX_INTO_CHAMBER + FUEL_INTO_CHAMBER)
    Chamber.mass_flow_out = Engine["mass_flow"]
    THRUST = Engine["thrust"]
    TOTAL_FLOW = OX_INTO_CHAMBER + FUEL_INTO_CHAMBER
    # FullFlow's SpecificImpulse; the flow floor only avoids dividing by zero
    # before ignition (the export reports Isp as 0 while the chamber is unlit).
    Isp = SpecificImpulse("Main Chamber Specific Impulse", Hopper, thrust=THRUST,
                          mass_flow=State.maximum(TOTAL_FLOW, 1.0e-6))


    # ---- Tracked States ---- #
    # Saved under "Hopper_Prop_System/transient/runs/base/tracks". The per-element
    # and per-node tracks are written by the builders above.
    Hopper.track("COPV Pressure [psia]", COPV.pressure / psi)
    Hopper.track("COPV Temperature [K]", COPVGas.temperature)
    Hopper.track("COPV Gas Mass [kg]", COPV_GAS_MASS)
    Hopper.track("COPV Density [kg/m3]", COPVGas.density)
    Hopper.track("COPV Enthalpy [J/kg]", COPVGas.enthalpy)

    for label, tank in [("LOX", LOX_TANK), ("Fuel", FUEL_TANK)]:
        parts = tanks[label]
        Hopper.track(f"{label} Tank Pressure [psia]", parts["ullage_gas"].pressure / psi)
        Hopper.track(f"{label} Liquid Mass [kg]", parts["liquid_node"].mass)
        Hopper.track(f"{label} Liquid Volume [L]", parts["liquid_node"].volume / liter)
        Hopper.track(f"{label} Liquid Density [kg/m3]", parts["liquid"].density)
        Hopper.track(f"{label} Liquid Enthalpy [J/kg]", parts["liquid"].enthalpy)
        Hopper.track(f"{label} Liquid Quality [-]", parts["liquid"].quality)
        Hopper.track(f"{label} Ullage Gas Mass [kg]", parts["ullage"].mass)
        Hopper.track(f"{label} Ullage Temperature [K]", parts["ullage_gas"].temperature)
        Hopper.track(f"{label} Ullage Density [kg/m3]", parts["ullage_gas"].density)
        Hopper.track(f"{label} Ullage Enthalpy [J/kg]", parts["ullage_gas"].enthalpy)
        Hopper.track(f"{label} Ullage Volume [L]", parts["ullage"].volume / liter)
        Hopper.track(f"{label} Ullage Inflow [kg/s]", parts["ullage"].mass_flow_in)
        Hopper.track(f"{label} Bang-Bang State", legs[label]["command"])
        Hopper.track(f"{label} Pressurant Mass Flow [kg/s]", legs[label]["mass_flow"])
        Hopper.track(f"{label} Tank Outflow [kg/s]", feeds[label]["tank_outflow"])
        for node_name, _, gas in legs[label]["nodes"]:
            Hopper.track(f"{node_name} Temperature [K]", gas.temperature)
            Hopper.track(f"{node_name} Density [kg/m3]", gas.density)

    Hopper.track("Mixture Ratio", MR)
    Hopper.track("Chamber Pressure [psia]", PC / psi)
    Hopper.track("Chamber Temperature [K]", Engine["gas"].temperature)
    Hopper.track("C* Efficiency", Engine["cstar_efficiency"])
    Hopper.track("Cf Efficiency", Engine["cf_efficiency"])
    Hopper.track("Nozzle Mass Flow [kg/s]", Chamber.mass_flow_out)
    Hopper.track("Ox Mass Flow [kg/s]", OX_INTO_CHAMBER)
    Hopper.track("Fuel Mass Flow [kg/s]", FUEL_INTO_CHAMBER)
    Hopper.track("Total Mass Flow [kg/s]", TOTAL_FLOW)
    Hopper.track("Thrust [N]", THRUST)
    Hopper.track("Target Thrust [N]", THRUST_TARGET)
    Hopper.track("Isp [s]", Isp.specific_impulse)
    # Delivered c* = Pc At / mdot (the floor only avoids dividing by zero before ignition).
    Hopper.track("C* [m/s]", PC * THROAT_AREA / State.maximum(Chamber.mass_flow_out, 1.0e-6))


    # ---- Solver ---- #
    # Implicit (backward Euler) transient. Each step solves all unknowns
    # together with SciPy least squares; a step whose residual exceeds rtol is
    # retried with half the time step.
    # Resolve physical travel with >=10 steps/full stroke and controller samples
    # with >=2 steps. This is resolution control, not a temporal error estimate.
    from Backend.Hopper_Config import integration_step
    resolved_dt = integration_step(globals(), dt)
    print(f"Network ready in {time.perf_counter() - build_start:.1f} s; "
          f"starting transient: {t_final:g} simulated seconds at <= {resolved_dt:g} s per step.", flush=True)
    if resolved_dt < dt:
        print(f"Valve/controller resolution: integration dt capped at {resolved_dt:g} s (requested {dt:g} s).")
    from Backend.Hopper_Solver import transient_solver
    transient_solver(Hopper, max_evaluations=PROP_MAX_NFEV).solve(
        dt=resolved_dt,
        t_final=t_final,
        save_dt=save_dt,
        filename=filename,
        rtol=rtol,
        jacobian_method=jacobian_method,
        verbose=verbose,
        statistics=statistics,
    )
    from Backend.Hopper_Report import write_run_metadata
    write_run_metadata(filename, targets, resolved_dt)

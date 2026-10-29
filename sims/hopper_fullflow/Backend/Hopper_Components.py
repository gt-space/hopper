"""Custom FullFlow components and the shared engine-performance block.
"""

import math

from fullflow import Component, IsentropicNozzle, Map, Sequence, State
from scipy.constants import g, psi

from HOPPER_DATA import *


class LinearizedNozzle(IsentropicNozzle):
    """Used to make the runtime acceptable for unchoked flow at the band edge (Pc = Pa + 1 psi).
    """

    def evaluate_states(self):
        # Full parent calculation first: exit state, choking, mass flow.
        super().evaluate_states()

        dp = self.upstream_total_pressure.value - self.ambient_pressure.value
        if dp >= NOZZLE_LINEAR_PRESSURE_DROP:
            return

        # Unchoked nozzle flow at the band edge (Pc = ambient + band), using the
        # parent's own flow-parameter function, then scaled linearly to dP.
        gamma = self.specific_heat_ratio.value
        rt = self.gas_constant.value * self.upstream_total_temperature.value
        edge_pressure = self.ambient_pressure.value + NOZZLE_LINEAR_PRESSURE_DROP
        flow_parameter = self._flow_parameter_from_pressure_ratio(edge_pressure / self.ambient_pressure.value, gamma)
        edge_mass_flow = (
            flow_parameter * self.expansion_ratio.value * edge_pressure * self.throat_area.value / math.sqrt(rt)
        )
        self.mass_flow.value = edge_mass_flow * dp / NOZZLE_LINEAR_PRESSURE_DROP


class DrainingTank(Component):
    """Liquid tank emptying through a bottom outlet (LOX drain calculator).

    State: liquid_mass, d/dt = -mass_flow_out. The outlet sees the ullage
    pressure plus the liquid head, rho g (h + outlet_drop), with the liquid
    height h = mass / (density * cross_sectional_area).
    """

    def __init__(self, name, network, liquid_mass: State, ullage_pressure, density, cross_sectional_area,
                 outlet_drop=0.0, mass_flow_out=0.0, outlet_pressure: State | None = None,
                 liquid_height: State | None = None):
        self.setup()

    def evaluate_states(self):
        height = max(self.liquid_mass.value, 0.0) / (self.density.value * self.cross_sectional_area.value)
        self.liquid_height.value = height
        self.outlet_pressure.value = (self.ullage_pressure.value
                                      + self.density.value * g * (height + self.outlet_drop.value))

    @property
    def dynamics(self):
        return [(self.liquid_mass, -self.mass_flow_out.value)]


def _smooth_unit_clip(x, width=0.02):
    """clip(x, 0, 1) with rounded corners: each corner is a sqrt-softplus of
    half-width ~width, within width/2 of the hard clip and equal to x inside."""
    def softplus(y):
        return 0.5 * (y + (y * y + width * width) ** 0.5)
    return softplus(x) - softplus(x - 1.0)


def efficiency(pc_psia, low, high, pc_low, pc_span):
    """eta = low + (high - low) * clip((Pc_psia - pc_low) / pc_span, 0, 1), as a plain number."""
    return low + (high - low) * min(max((pc_psia - pc_low) / pc_span, 0.0), 1.0)


def EnginePerformance(name, network, chamber_pressure, mixture_ratio,
                      chamber_map_filename=CHAMBER_MAP_FILENAME,
                      ambient_pressure=AMBIENT_PRESSURE,
                      sampled_efficiency=False):
    """Chamber map, nozzle, efficiencies, real flow and thrust.

    sampled_efficiency=False: the efficiencies are expressions of Pc (steady
    solves: thrust targeter, sweeps). True (the transient): each efficiency
    is a FullFlow function Sequence evaluated at the previous accepted step's
    Pc, like a sampled controller. Making the transient's residuals depend on
    Pc through the efficiency curve stalls FullFlow's Newton path at ignition
    (every step falls back to SciPy TRF); a one-step lag (0.05 s) changes eta
    by well under 0.1 %.
    """
    # Equilibrium chamber products: temperature, gamma, gas constant, density.
    gas = Map.from_hdf5(
        f"{name} Gas",
        network,
        chamber_map_filename,
        group="chamber",
        inputs={
            "chamber_pressure": chamber_pressure,
            "mixture_ratio": mixture_ratio,
        },
        extrapolate=True,   # tolerate solver trial points just outside the map
    )

    # Efficiency ramps linearly from LOW at PC_LOW to HIGH at PC_LOW + SPAN, flat outside:
    #   eta = LOW + (HIGH - LOW) * clip((Pc_psia - PC_LOW) / SPAN, 0, 1)
    # (vehicle_sim.py used max(1, ...), which pinned eta at HIGH below the top of the ramp.)
    # In the steady targeter the clip's corners are rounded (_smooth_unit_clip) so
    # the Pc solve has a continuous slope.
    if sampled_efficiency:
        cstar_efficiency, cf_efficiency = State(CSTAR_EFFICIENCY_LOW), State(CF_EFFICIENCY_LOW)
        Sequence(f"{name} C* Efficiency", network, target=cstar_efficiency, inputs=[chamber_pressure],
                 function=lambda t, pc: efficiency(pc / psi, CSTAR_EFFICIENCY_LOW, CSTAR_EFFICIENCY_HIGH,
                                                   CSTAR_EFFICIENCY_PC_LOW, CSTAR_EFFICIENCY_PC_SPAN))
        Sequence(f"{name} Cf Efficiency", network, target=cf_efficiency, inputs=[chamber_pressure],
                 function=lambda t, pc: efficiency(pc / psi, CF_EFFICIENCY_LOW, CF_EFFICIENCY_HIGH,
                                                   CF_EFFICIENCY_PC_LOW, CF_EFFICIENCY_PC_SPAN))
        return _engine_outputs(name, network, chamber_pressure, gas, ambient_pressure,
                               cstar_efficiency, cf_efficiency)

    pc_psia = chamber_pressure / psi
    cstar_efficiency = CSTAR_EFFICIENCY_LOW + (CSTAR_EFFICIENCY_HIGH - CSTAR_EFFICIENCY_LOW) * _smooth_unit_clip(
        (pc_psia - CSTAR_EFFICIENCY_PC_LOW) / CSTAR_EFFICIENCY_PC_SPAN)
    cf_efficiency = CF_EFFICIENCY_LOW + (CF_EFFICIENCY_HIGH - CF_EFFICIENCY_LOW) * _smooth_unit_clip(
        (pc_psia - CF_EFFICIENCY_PC_LOW) / CF_EFFICIENCY_PC_SPAN)
    return _engine_outputs(name, network, chamber_pressure, gas, ambient_pressure,
                           cstar_efficiency, cf_efficiency)


def _engine_outputs(name, network, chamber_pressure, gas, ambient_pressure, cstar_efficiency, cf_efficiency):
    # Ideal (100 % c*) nozzle: flow, exit velocity and exit pressure at the
    # ambient back pressure (sea-level test conditions, so Pe < Pa is possible).
    nozzle = LinearizedNozzle(
        f"{name} Nozzle",
        network,
        upstream_total_pressure=chamber_pressure,
        upstream_total_temperature=gas.temperature,
        ambient_pressure=ambient_pressure,
        specific_heat_ratio=gas.gamma,
        gas_constant=gas.gas_constant,
        throat_area=THROAT_AREA,
        expansion_ratio=EXPANSION_RATIO,
    )

    return {
        "gas": gas,
        "nozzle": nozzle,
        "cstar_efficiency": cstar_efficiency,
        "cf_efficiency": cf_efficiency,
        # Real chamber outflow: more propellant than the ideal nozzle passes,
        # because the real c* is only eta_c* of ideal.
        "mass_flow": nozzle.mass_flow / cstar_efficiency,
        # Momentum + pressure thrust of the ideal nozzle, scaled by eta_cf.
        "thrust": cf_efficiency * (
            nozzle.mass_flow * nozzle.exit_velocity
            + (nozzle.exit_static_pressure - ambient_pressure) * EXIT_AREA
        ),
    }

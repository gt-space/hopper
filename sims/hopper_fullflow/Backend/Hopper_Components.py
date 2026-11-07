"""Custom FullFlow components and the shared engine-performance block.
"""

import math

from fullflow import Component, IsentropicNozzle, Map, Sequence, State
from scipy.constants import g

from HOPPER_DATA import *
from Backend.Hopper_Inputs import efficiency_value
from Backend.Hopper_Control import FunctionLookup


def engine_efficiencies(chamber_pressure, mixture_ratio):
    return {"cstar": CSTAR_EFFICIENCY_SCALE * efficiency_value(CSTAR_EFFICIENCY_FUNCTION, CSTAR_EFFICIENCY_CSV,
                                                               chamber_pressure, mixture_ratio),
            "cf": CF_EFFICIENCY_SCALE * efficiency_value(CF_EFFICIENCY_FUNCTION, CF_EFFICIENCY_CSV,
                                                         chamber_pressure, mixture_ratio)}


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


def EnginePerformance(name, network, chamber_pressure, mixture_ratio,
                      chamber_map_filename=CHAMBER_MAP_FILENAME,
                      ambient_pressure=AMBIENT_PRESSURE,
                      sampled_efficiency=False):
    """Chamber map, nozzle, efficiencies, real flow and thrust.

    sampled_efficiency=False: the efficiencies follow Pc in the same solve (steady
    solves: thrust targeter, sweeps). True (the transient): each efficiency is a
    FullFlow function Sequence evaluated at the previous accepted step's Pc, like a
    sampled controller. Making the transient's residuals depend on Pc through the
    efficiency curve stalls FullFlow's Newton path at ignition (every step falls back
    to SciPy TRF); a one-step lag (0.05 s) changes eta by well under 0.1 %.
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

    if sampled_efficiency:
        start = engine_efficiencies(ambient_pressure, MIXTURE_RATIO)
        cstar_efficiency, cf_efficiency = State(start["cstar"]), State(start["cf"])
        for label, target, key in [("C*", cstar_efficiency, "cstar"), ("Cf", cf_efficiency, "cf")]:
            Sequence(f"{name} {label} Efficiency", network, target=target, inputs=[chamber_pressure, mixture_ratio],
                     function=lambda t, pc, mr, key=key: engine_efficiencies(pc, mr)[key])
    else:
        efficiencies = FunctionLookup(f"{name} Efficiencies", network, engine_efficiencies,
                              chamber_pressure=chamber_pressure, mixture_ratio=mixture_ratio,
                              outputs=("cstar", "cf"))
        cstar_efficiency, cf_efficiency = efficiencies.cstar, efficiencies.cf

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

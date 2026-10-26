"""Custom FullFlow components and the shared engine-performance block.
"""

import math

from fullflow import IsentropicNozzle, Map, State
from scipy.constants import psi

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


def EnginePerformance(name, network, chamber_pressure, mixture_ratio,
                      chamber_map_filename=CHAMBER_MAP_FILENAME,
                      ambient_pressure=AMBIENT_PRESSURE):
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

    # vehicle_sim.dynamic_eta_cstar / dynamic_eta_cf, written exactly as there:
    #   eta = LOW + (HIGH - LOW) * max(1, (Pc_psia - PC_LOW) / SPAN)
    # (see HOPPER_DATA for what max(1, ...) does to the curves).
    pc_psia = chamber_pressure / psi
    cstar_efficiency = CSTAR_EFFICIENCY_LOW + (CSTAR_EFFICIENCY_HIGH - CSTAR_EFFICIENCY_LOW) * State.maximum(
        1.0, (pc_psia - CSTAR_EFFICIENCY_PC_LOW) / CSTAR_EFFICIENCY_PC_SPAN
    )
    cf_efficiency = CF_EFFICIENCY_LOW + (CF_EFFICIENCY_HIGH - CF_EFFICIENCY_LOW) * State.maximum(
        1.0, (pc_psia - CF_EFFICIENCY_PC_LOW) / CF_EFFICIENCY_PC_SPAN
    )

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

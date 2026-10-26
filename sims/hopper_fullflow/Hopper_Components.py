"""Custom FullFlow components and the shared engine-performance block.

Why custom components are needed
--------------------------------
FullFlow's IsentropicNozzle uses the textbook ideal-gas flow law. Near zero
pressure difference, flow goes as sqrt(dP), whose slope is infinite at
dP = 0. Hopper's main chamber sits at exactly that point before start-up and
after shutdown, and the solver stalls there (it cannot form a usable
Jacobian). LinearizedNozzle below is the same FullFlow relation, with the flow
made linear in dP inside a small band around zero. The flow is continuous at
the band edge, and outside the band the result is identical to the stock
component.

Engine performance
------------------
EnginePerformance() builds the chamber gas, nozzle and thrust from (Pc, MR).
The thrust targeter and the prop system both call it, so the flows the targeter
asks for are exactly what the transient engine needs to make that thrust.
"""

import math

from fullflow import IsentropicNozzle, Map, State
from scipy.constants import psi

from HOPPER_DATA import *


class LinearizedNozzle(IsentropicNozzle):
    """FullFlow's IsentropicNozzle with its mass flow linearised near ambient.

    While the chamber is within NOZZLE_LINEAR_PRESSURE_DROP of ambient (before
    ignition and after tail-off), the nozzle flow is the band-edge flow scaled
    linearly with (Pc - ambient). Everything else (choking, exit Mach, exit
    pressure and velocity, shock checks) is the parent component unchanged.
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
    """Chamber gas, nozzle flow and thrust at (Pc, MR).

    chamber_pressure and mixture_ratio are States (or live State expressions);
    everything returned stays connected to them.

    Physics, matching vehicle_sim.Engine and its targeter:
      gas          equilibrium RP-1/LOX products at (Pc, MR) from the chamber map
      mdot_ideal   isentropic-nozzle flow for that gas (= Pc At / c*_ideal when choked)
      mdot         mdot_ideal / eta_c*                (c* efficiency, as in Vespula)
      thrust       eta_cf * (mdot_ideal * Ve + (Pe - Pa) * Ae)
                   = eta_cf * Cf_ambient * Pc * At    (vehicle_sim: Isp_amb * eta_c* * eta_cf * g * mdot)

    Returns a dict of the map, nozzle, both efficiencies, mdot and thrust.
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

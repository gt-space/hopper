"""Thrust-profile targeter: thrust [N] -> required Pc and propellant flows.

Replaces vehicle_sim.ThrustProfileTargeter. That class solved, for each point
of the thrust profile, the chamber pressure at which RocketCEA predicts the
target thrust at MR 1.8, and then split the resulting mass flow into ox and
fuel. This does the same thing, but each point is a steady FullFlow solve of
the same EnginePerformance block the transient uses, so the targets and the
transient engine are always consistent (no RocketCEA/FullFlow mismatch).
"""

from fullflow import Balance, Network, State, SteadyState
from scipy.constants import *

from HOPPER_DATA import *
from Hopper_Components import EnginePerformance


def HopperThrustTargets(thrust_profile: dict,
                        mixture_ratio: float = MIXTURE_RATIO,
                        chamber_map_filename: str = CHAMBER_MAP_FILENAME):
    """Required Pc and propellant flows for every point of a thrust profile.

    thrust_profile maps profile time [s] to thrust [N] (e.g. THRUST_CURVE_1DOF).

    Returns a list of dicts, one per profile point, with:
      time              profile time [s]
      thrust            target thrust [N]
      chamber_pressure  Pc that gives that thrust at the target MR [Pa]
      ox_mass_flow      required LOX flow [kg/s]
      fuel_mass_flow    required RP-1 flow [kg/s]
    """

    # A tiny steady network: one unknown (Pc) and one equation (thrust = target).
    Targeter = Network("Hopper Thrust Targeter")

    chamber_pressure = State(250 * psi)                          # initial guess; solved below
    target_thrust = State(next(iter(thrust_profile.values())))   # updated for each profile point

    Engine = EnginePerformance(
        "Targeter Chamber",
        Targeter,
        chamber_pressure=chamber_pressure,
        mixture_ratio=State(mixture_ratio),   # fixed at the target MR
        chamber_map_filename=chamber_map_filename,
    )

    # The solver varies Pc until thrust matches the target. The residual is the
    # thrust error divided by 100 N.
    Balance(
        "Thrust Target",
        Targeter,
        variable=chamber_pressure,
        function=(Engine["thrust"] - target_thrust) / 100.0,
    )

    targets = []
    for time, thrust in thrust_profile.items():
        target_thrust.value = thrust
        # rtol = 1e-4 on (error / 100 N) accepts a 0.01 N thrust miss. Tight
        # ftol/xtol/gtol keep SciPy iterating until FullFlow's rtol test passes.
        # Each solve starts from the previous point's Pc.
        SteadyState(Targeter).solve(rtol=1e-4, ftol=1e-12, xtol=1e-12, gtol=1e-12)

        # Total flow out of the chamber at that Pc, split at the target MR.
        mass_flow = Engine["mass_flow"].value
        fuel_mass_flow = mass_flow / (1.0 + mixture_ratio)
        targets.append({
            "time": time,
            "thrust": thrust,
            "chamber_pressure": chamber_pressure.value,
            "ox_mass_flow": mass_flow - fuel_mass_flow,
            "fuel_mass_flow": fuel_mass_flow,
        })

    return targets

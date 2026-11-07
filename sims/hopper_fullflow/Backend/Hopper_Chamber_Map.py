"""Combustion-product property map for the main chamber.
"""

from thermoprop import Propellant, Reactants, Equilibrium
import fullplot as fplt
from scipy.constants import *

from HOPPER_DATA import *


def HopperChamberMap(filename: str = CHAMBER_MAP_FILENAME,
                     lox_temperature: float = LOX_TANK.liquid_temperature,
                     fuel_temperature: float = FUEL_TANK.liquid_temperature,
                     chamber_pressure_lower_limit: float = CHAMBER_MAP_PRESSURE_LOWER_LIMIT,
                     chamber_pressure_upper_limit: float = CHAMBER_MAP_PRESSURE_UPPER_LIMIT,
                     chamber_pressure_data_points: int = CHAMBER_MAP_PRESSURE_POINTS,
                     mixture_ratio_lower_limit: float = CHAMBER_MAP_MIXTURE_RATIO_LOWER_LIMIT,
                     mixture_ratio_upper_limit: float = CHAMBER_MAP_MIXTURE_RATIO_UPPER_LIMIT,
                     mixture_ratio_data_points: int = CHAMBER_MAP_MIXTURE_RATIO_POINTS,
                     fuel_name: str = FUEL,
                     oxidizer_name: str = OXIDIZER,
                     ):
    """Build the HDF5 map (default grid 31 x 36 equilibrium solves, ~2 minutes).

    The grid and propellants come from HOPPER_DATA.py. The map is stored in
    group "chamber" of `filename`.
    """

    # Reactants at their tank temperatures (enthalpy reference for the
    # adiabatic flame temperature).
    fuel = Propellant(
        fuel_name,
        temperature=fuel_temperature,
    )

    oxidizer = Propellant(
        oxidizer_name,
        temperature=lox_temperature,
    )


    def chamber_products(chamber_pressure, mixture_ratio):
        """Adiabatic equilibrium products at one grid point."""

        reactants = Reactants(
            fuels=fuel,
            oxidizers=oxidizer,
            mixture_ratio=mixture_ratio,
        )

        gas = Equilibrium(
            reactants=reactants,
            pressure=chamber_pressure,
        )

        # Exactly what the chamber Volume and the nozzle need.
        return {
            "temperature": gas.temperature,
            "gamma": gas.gamma,
            "gas_constant": gas.gas_constant,
            "density": gas.density,
        }


    fplt.generate_map(
        filename,
        group="chamber",
        axes=[
            fplt.Axis.linear(
                "chamber_pressure",
                start=chamber_pressure_lower_limit,
                stop=chamber_pressure_upper_limit,
                count=chamber_pressure_data_points,
                units="Pa",
            ),
            fplt.Axis.linear(
                "mixture_ratio",
                start=mixture_ratio_lower_limit,
                stop=mixture_ratio_upper_limit,
                count=mixture_ratio_data_points,
            ),
        ],
        evaluate=chamber_products,
        overwrite=True,
        raise_errors=True,   # stop if any grid point fails, rather than leave a hole
    )

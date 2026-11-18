"""Preliminary steady-state ORV sizing for oxidizer B/B valves failed open.

The model screens gas-only relief of GN2 inflow plus provisional simultaneous
GOX generation. It is not a transient model or certified relief-capacity
method. The configured multi-restriction path uses the repository SeriesCdA
convention, and ideal-gas isentropic flow selects choked or unchoked behavior.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

from ORV_config import DEFAULT_HOPPER_CONFIG, ORVConfigError, load_orv_config

PSI_TO_PA = 6894.757293168
INCH_TO_M = 0.0254
MM2_TO_M2 = 1e-6
M2_TO_MM2 = 1e6


def ideal_gas_properties(fluid: str, temperature_k: float) -> dict[str, float]:
    """Calculate ideal-gas R, cp, cv, and gamma from CoolProp fluid data."""
    if PropsSI is None:
        raise RuntimeError("CoolProp is required to calculate gas properties")
    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("Gas temperature must be finite and positive")
    gas_constant = float(PropsSI("GAS_CONSTANT", fluid) / PropsSI("MOLARMASS", fluid))
    cp = float(PropsSI("CP0MASS", "T", temperature_k, "P", 101325.0, fluid))
    cv = cp - gas_constant
    return {
        "specific_gas_constant_j_kg_k": gas_constant,
        "cp_j_kg_k": cp,
        "cv_j_kg_k": cv,
        "gamma": cp / cv,
        "temperature_k": temperature_k,
    }


def ideal_gas_mixture_properties(
    nitrogen_mass_flow_kg_s: float,
    oxygen_mass_flow_kg_s: float,
    temperature_k: float,
) -> dict[str, float]:
    """Mass-weight component ideal-gas properties at the relieving temperature."""
    flows = np.asarray([nitrogen_mass_flow_kg_s, oxygen_mass_flow_kg_s], dtype=float)
    if np.any(~np.isfinite(flows)) or np.any(flows < 0.0) or float(np.sum(flows)) <= 0.0:
        raise ValueError("Mixture mass flows must be finite, nonnegative, and total positive")
    fractions = flows / np.sum(flows)
    nitrogen = ideal_gas_properties("Nitrogen", temperature_k)
    oxygen = ideal_gas_properties("Oxygen", temperature_k)
    gas_constant = float(fractions[0] * nitrogen["specific_gas_constant_j_kg_k"] + fractions[1] * oxygen["specific_gas_constant_j_kg_k"])
    cp = float(fractions[0] * nitrogen["cp_j_kg_k"] + fractions[1] * oxygen["cp_j_kg_k"])
    cv = cp - gas_constant
    return {
        "nitrogen_mass_fraction": float(fractions[0]),
        "oxygen_mass_fraction": float(fractions[1]),
        "specific_gas_constant_j_kg_k": gas_constant,
        "cp_j_kg_k": cp,
        "cv_j_kg_k": cv,
        "gamma": cp / cv,
        "temperature_k": temperature_k,
    }


def oxygen_saturation_properties(pressure_pa: float, fluid: str = "Oxygen") -> dict[str, float]:
    """Return saturated-oxygen properties below the critical pressure."""
    if PropsSI is None:
        raise RuntimeError("CoolProp is required to calculate oxygen properties")
    critical_pressure = float(PropsSI("PCRIT", fluid))
    if pressure_pa >= critical_pressure:
        raise ValueError("Oxygen saturation properties are undefined at/above critical pressure")
    liquid_h = float(PropsSI("H", "P", pressure_pa, "Q", 0, fluid))
    vapor_h = float(PropsSI("H", "P", pressure_pa, "Q", 1, fluid))
    return {
        "saturation_temperature_k": float(PropsSI("T", "P", pressure_pa, "Q", 0, fluid)),
        "latent_heat_j_kg": vapor_h - liquid_h,
    }


def determine_gox_generation(config: dict[str, Any], pressure_psia: float) -> dict[str, Any]:
    """Use Qdot/h_fg only for saturated LOX; otherwise retain provisional GOX."""
    gox = config["gox"]
    provisional = gox["provisional_generation_kg_s"]
    try:
        saturation = oxygen_saturation_properties(pressure_psia * PSI_TO_PA, gox["fluid_name"])
    except (RuntimeError, ValueError) as exc:
        return {"mass_flow_kg_s": provisional, "source": "provisional_unverified", "status": f"provisional GOX retained; saturation check unavailable: {exc}", "saturation_temperature_k": None, "subcooling_k": None}
    subcooling = saturation["saturation_temperature_k"] - config["tank"]["liquid_temperature_k"]
    heat_rate = gox["effective_heat_transfer_rate_w"]
    if heat_rate is not None and subcooling <= 0.1:
        mass_flow = heat_rate * gox["heat_to_vaporization_fraction"] / saturation["latent_heat_j_kg"]
        source = "calculated_from_heat"
        status = "Qdot_effective/h_fg for approximately saturated LOX"
    else:
        mass_flow = provisional
        source = "provisional_unverified"
        if heat_rate is None:
            status = "provisional GOX retained; no validated heat input"
        else:
            status = "provisional GOX retained; LOX is subcooled and sensible-heating inventory/time are unavailable"
    return {"mass_flow_kg_s": mass_flow, "source": source, "status": status, "saturation_temperature_k": saturation["saturation_temperature_k"], "subcooling_k": subcooling, "latent_heat_j_kg": saturation["latent_heat_j_kg"]}


def ideal_gas_enthalpy_change(fluid: str, inlet_temperature_k: float, outlet_temperature_k: float) -> float:
    """Return h(out)-h(in) [J/kg] by integrating CoolProp ideal-gas cp0."""
    for value, label in ((inlet_temperature_k, "inlet temperature"), (outlet_temperature_k, "outlet temperature")):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{label} must be in kelvin and greater than zero")
    if PropsSI is None:
        raise RuntimeError("CoolProp is required to calculate ideal-gas enthalpy")
    if inlet_temperature_k == outlet_temperature_k:
        return 0.0
    temperatures = np.linspace(inlet_temperature_k, outlet_temperature_k, 129)
    heat_capacities = np.asarray([
        PropsSI("CP0MASS", "T", temperature, "P", 101325.0, fluid)
        for temperature in temperatures
    ], dtype=float)
    return float(np.trapezoid(heat_capacities, temperatures))


def solve_well_mixed_ullage_temperature(
    gn2_mass_flow_kg_s: float,
    gn2_inlet_temperature_k: float,
    gox_mass_flow_kg_s: float,
    gox_inlet_temperature_k: float,
    liquid_temperature_k: float,
    liquid_to_ullage_ua_w_k: float = 0.0,
    wall_temperature_k: float | None = None,
    wall_to_ullage_ua_w_k: float = 0.0,
    net_external_sensible_heat_rate_w: float = 0.0,
) -> dict[str, float]:
    """Solve a steady, well-mixed ullage enthalpy balance.

    GOX enters as already-generated vapor, so latent heat is not added here.
    Positive heat rates are into the gas. The ORV outlet carries both species at
    the solved mixed temperature and their physical mass flows (no design factor).
    """
    for value, label in (
        (gn2_mass_flow_kg_s, "GN2 mass flow"), (gox_mass_flow_kg_s, "GOX mass flow"),
        (liquid_to_ullage_ua_w_k, "liquid UA"), (wall_to_ullage_ua_w_k, "wall UA"),
    ):
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{label} must be finite and nonnegative")
    if gn2_mass_flow_kg_s + gox_mass_flow_kg_s <= 0.0:
        raise ValueError("Total ullage gas flow must be positive")
    for value, label in ((gn2_inlet_temperature_k, "GN2 inlet temperature"), (gox_inlet_temperature_k, "GOX inlet temperature"), (liquid_temperature_k, "liquid temperature")):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{label} must be in kelvin and greater than zero")
    if wall_to_ullage_ua_w_k > 0.0 and wall_temperature_k is None:
        raise ValueError("Wall temperature is required when wall UA is positive")
    if wall_temperature_k is not None and (not math.isfinite(wall_temperature_k) or wall_temperature_k <= 0.0):
        raise ValueError("Wall temperature must be in kelvin and greater than zero")
    if not math.isfinite(net_external_sensible_heat_rate_w):
        raise ValueError("External sensible heat rate must be finite")

    def heat_into_gas(temperature_k: float) -> float:
        liquid_heat = liquid_to_ullage_ua_w_k * (liquid_temperature_k - temperature_k)
        wall_heat = 0.0 if wall_temperature_k is None else wall_to_ullage_ua_w_k * (wall_temperature_k - temperature_k)
        return liquid_heat + wall_heat + net_external_sensible_heat_rate_w

    def residual(temperature_k: float) -> float:
        outlet_minus_inlet = (
            gn2_mass_flow_kg_s * ideal_gas_enthalpy_change("Nitrogen", gn2_inlet_temperature_k, temperature_k)
            + gox_mass_flow_kg_s * ideal_gas_enthalpy_change("Oxygen", gox_inlet_temperature_k, temperature_k)
        )
        return outlet_minus_inlet - heat_into_gas(temperature_k)

    boundary_temperatures = [gn2_inlet_temperature_k, gox_inlet_temperature_k, liquid_temperature_k]
    if wall_temperature_k is not None:
        boundary_temperatures.append(wall_temperature_k)
    lower = max(1.0, min(boundary_temperatures) * 0.5)
    upper = max(boundary_temperatures) * 1.5
    lower_residual, upper_residual = residual(lower), residual(upper)
    while lower_residual * upper_residual > 0.0 and upper < 5000.0:
        lower = max(1.0, lower * 0.5)
        upper *= 1.5
        lower_residual, upper_residual = residual(lower), residual(upper)
    if lower_residual * upper_residual > 0.0:
        raise RuntimeError("Could not bracket the ullage energy-balance temperature")
    for _ in range(100):
        temperature = 0.5 * (lower + upper)
        middle_residual = residual(temperature)
        if abs(middle_residual) <= 1e-7 or upper - lower <= 1e-8:
            break
        if lower_residual * middle_residual <= 0.0:
            upper, upper_residual = temperature, middle_residual
        else:
            lower, lower_residual = temperature, middle_residual
    else:
        raise RuntimeError("Ullage energy-balance solver did not converge")
    return {
        "temperature_k": temperature,
        "net_heat_transfer_w": heat_into_gas(temperature),
        "energy_residual_w": residual(temperature),
    }


def calculate_orv_relieving_temperature(
    config: dict[str, Any],
    gn2_mass_flow_kg_s: float,
    gox: dict[str, Any],
) -> dict[str, Any]:
    """Select calculated temperature only when a thermal closure is available."""
    thermal = config["thermal"]
    gox_inlet_temperature = gox.get("saturation_temperature_k")
    adiabatic_limit = None
    if gox_inlet_temperature is not None:
        adiabatic_limit = solve_well_mixed_ullage_temperature(
            gn2_mass_flow_kg_s, min(config["source"]["temperature_k"]),
            gox["mass_flow_kg_s"], gox_inlet_temperature,
            config["tank"]["liquid_temperature_k"],
        )

    has_ua_closure = (
        thermal["liquid_to_ullage_ua_w_k"] is not None
        and thermal["wall_to_ullage_ua_w_k"] is not None
        and thermal["tank_wall_temperature_k"] is not None
    )
    has_direct_heat_closure = thermal["net_external_sensible_heat_rate_w"] is not None
    missing = []
    if not has_ua_closure and not has_direct_heat_closure:
        missing.extend((
            "liquid-to-ullage UA or validated equivalent heat-transfer model",
            "wall-to-ullage UA and relief-event wall temperature",
            "net sensible heat rate to ullage, distinct from vaporization energy",
        ))
    if gox_inlet_temperature is None:
        missing.append("GOX inlet thermodynamic state")

    calculated = None
    if not missing:
        calculated = solve_well_mixed_ullage_temperature(
            gn2_mass_flow_kg_s, min(config["source"]["temperature_k"]),
            gox["mass_flow_kg_s"], gox_inlet_temperature,
            config["tank"]["liquid_temperature_k"],
            liquid_to_ullage_ua_w_k=thermal["liquid_to_ullage_ua_w_k"] or 0.0,
            wall_temperature_k=thermal["tank_wall_temperature_k"],
            wall_to_ullage_ua_w_k=thermal["wall_to_ullage_ua_w_k"] or 0.0,
            net_external_sensible_heat_rate_w=thermal["net_external_sensible_heat_rate_w"] or 0.0,
        )
        selected_temperature = calculated["temperature_k"]
        method = "calculated_steady_well_mixed_enthalpy_balance"
        status = "calculated from physical inflow enthalpies and supplied thermal closure"
    else:
        selected_temperature = thermal["provisional_relieving_temperature_k"]
        method = "provisional_assumption_missing_thermal_closure"
        status = "not calculated; 150 K configured ullage assumption retained"
    implied_selected_heat = None
    if gox_inlet_temperature is not None:
        implied_selected_heat = (
            gn2_mass_flow_kg_s * ideal_gas_enthalpy_change("Nitrogen", min(config["source"]["temperature_k"]), selected_temperature)
            + gox["mass_flow_kg_s"] * ideal_gas_enthalpy_change("Oxygen", gox_inlet_temperature, selected_temperature)
        )
    return {
        "selected_temperature_k": selected_temperature,
        "method": method,
        "status": status,
        "calculated": calculated,
        "adiabatic_limiting_case": adiabatic_limit,
        "net_heat_transfer_w": None if calculated is None else calculated["net_heat_transfer_w"],
        "implied_heat_transfer_at_selected_temperature_w": implied_selected_heat,
        "missing_inputs": missing,
        "sensitivity_min_temperature_k": config["tank"]["liquid_temperature_k"],
        "sensitivity_max_temperature_k": min(config["source"]["temperature_k"]),
        "gox_inlet_temperature_k": gox_inlet_temperature,
    }


def series_equivalent_cda(areas_m2: Iterable[float]) -> float:
    """Repository FullFlow SeriesCdA convention for restrictions in series."""
    areas = np.asarray(tuple(areas_m2), dtype=float)
    if areas.ndim != 1 or areas.size == 0 or np.any(~np.isfinite(areas)) or np.any(areas <= 0.0):
        raise ValueError("Every effective CdA must be finite and positive")
    return float(1.0 / np.sqrt(np.sum(1.0 / areas**2)))


def compressible_mass_flow(
    upstream_pressure_pa: Any,
    upstream_temperature_k: Any,
    downstream_pressure_pa: Any,
    effective_cda_m2: Any,
    gas_constant_j_kg_k: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Ideal-gas isentropic orifice flow using absolute SI inputs."""
    upstream, temperature, downstream, cda = np.broadcast_arrays(
        np.asarray(upstream_pressure_pa, dtype=float), np.asarray(upstream_temperature_k, dtype=float),
        np.asarray(downstream_pressure_pa, dtype=float), np.asarray(effective_cda_m2, dtype=float),
    )
    for values, label in ((upstream, "upstream pressure"), (downstream, "downstream pressure"), (temperature, "temperature"), (cda, "CdA")):
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError(f"{label} must be finite and positive")
    if not math.isfinite(gas_constant_j_kg_k) or gas_constant_j_kg_k <= 0.0:
        raise ValueError("Specific gas constant must be finite and positive")
    if not math.isfinite(gamma) or gamma <= 1.0:
        raise ValueError("Specific-heat ratio must be finite and greater than one")

    ratio = downstream / upstream
    critical_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    forward = upstream > downstream
    choked = forward & (ratio <= critical_ratio)
    unchoked = forward & ~choked
    mass_flow = np.zeros(upstream.shape, dtype=float)
    choked_factor = np.sqrt(gamma / (gas_constant_j_kg_k * temperature)) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    mass_flow[choked] = cda[choked] * upstream[choked] * choked_factor[choked]
    term = ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma)
    unchoked_factor = np.sqrt(np.maximum(0.0, 2.0 * gamma / (gas_constant_j_kg_k * temperature * (gamma - 1.0)) * term))
    mass_flow[unchoked] = cda[unchoked] * upstream[unchoked] * unchoked_factor[unchoked]
    return mass_flow, choked


def build_failed_open_path(config: dict[str, Any]) -> dict[str, float]:
    """Build the configured COPV-to-LOX-ullage steady effective restriction."""
    line = config["failure"]["pressurant_line"]
    diameter_in = line.get("inner_diameter_in")
    if diameter_in is None:
        diameter_in = line["outer_diameter_in"] - 2.0 * line["wall_thickness_in"]
    diameter_m = diameter_in * INCH_TO_M
    length_m = line["length_in"] * INCH_TO_M
    roughness_m = line["roughness_in"] * INCH_TO_M
    area_m2 = math.pi * diameter_m**2 / 4.0
    friction = (-2.0 * math.log10(roughness_m / (3.7 * diameter_m))) ** -2.0
    line_cda_m2 = area_m2 / math.sqrt(1.0 + friction * length_m / diameter_m)
    parallel_bottle_outlets_m2 = sum(config["source"]["outlet_cda_mm2"]) * MM2_TO_M2
    path_cda_m2 = series_equivalent_cda((parallel_bottle_outlets_m2, line_cda_m2, config["failure"]["obb_cda_mm2"] * MM2_TO_M2, config["failure"]["orifice_cda_mm2"] * MM2_TO_M2, config["failure"]["check_valve_cda_mm2"] * MM2_TO_M2))
    return {"line_inside_diameter_in": diameter_in, "line_friction_factor": friction, "line_cda_m2": line_cda_m2, "path_cda_m2": path_cda_m2}


def calculate_sizing(config: dict[str, Any]) -> dict[str, Any]:
    """Calculate the conservative MEOP-inflow / exploratory-relief-pressure case."""
    if not config["sizing_ready"]:
        raise ORVConfigError("Numerical ORV sizing is blocked")
    path = build_failed_open_path(config)
    source_pressure_pa = max(config["source"]["pressure_psia"]) * PSI_TO_PA
    source_temperature_k = min(config["source"]["temperature_k"])
    gn2 = ideal_gas_properties(config["pressurant"]["name"], source_temperature_k)
    relief_pressure_psia = config["orv"]["evaluated_relieving_pressure_psia"]
    pressures_psia = np.linspace(config["tank"]["meop_psia"], relief_pressure_psia, config["analysis"]["pressure_sweep_points"])
    inflow, inflow_choked = compressible_mass_flow(source_pressure_pa, source_temperature_k, pressures_psia * PSI_TO_PA, path["path_cda_m2"], gn2["specific_gas_constant_j_kg_k"], gn2["gamma"])
    maximum_index = int(np.argmax(inflow))
    analytical_maximum_inflow = float(inflow[maximum_index])
    maximum_inflow = config["failure"]["simulation_maximum_gn2_inflow_kg_s"]
    analytical_difference_percent = 100.0 * (
        analytical_maximum_inflow - maximum_inflow
    ) / maximum_inflow
    maximum_inflow_pressure = float(pressures_psia[maximum_index])

    gox = determine_gox_generation(config, relief_pressure_psia)
    factor = config["failure"]["gn2_mass_flow_safety_factor"]
    design_gn2 = factor * maximum_inflow
    combined_design_flow = design_gn2 + gox["mass_flow_kg_s"]
    thermal = calculate_orv_relieving_temperature(config, maximum_inflow, gox)
    relieving_temperature_k = thermal["selected_temperature_k"]
    # The safety factor is a design allowance, not a change in physical mixture
    # composition; use the raw simultaneous GN2/GOX loads for mixture properties.
    mixture = ideal_gas_mixture_properties(maximum_inflow, gox["mass_flow_kg_s"], relieving_temperature_k)
    relieving_gn2 = ideal_gas_properties(config["pressurant"]["name"], relieving_temperature_k)
    gn2_unit_capacity, _ = compressible_mass_flow(relief_pressure_psia * PSI_TO_PA, relieving_temperature_k, config["ambient_pressure_psia"] * PSI_TO_PA, 1.0, relieving_gn2["specific_gas_constant_j_kg_k"], relieving_gn2["gamma"])
    required_gn2_cda = design_gn2 / float(gn2_unit_capacity) * M2_TO_MM2
    unit_capacity, relief_choked = compressible_mass_flow(relief_pressure_psia * PSI_TO_PA, relieving_temperature_k, config["ambient_pressure_psia"] * PSI_TO_PA, 1.0, mixture["specific_gas_constant_j_kg_k"], mixture["gamma"])
    required_cda = combined_design_flow / float(unit_capacity) * M2_TO_MM2
    required_diameter = math.sqrt(4.0 * required_cda / math.pi)
    candidate_cda = config["orv"]["candidate_cda_mm2"]
    candidate_capacity = float(unit_capacity) * candidate_cda * MM2_TO_M2

    sweep_unit_capacity, sweep_relief_choked = compressible_mass_flow(pressures_psia * PSI_TO_PA, relieving_temperature_k, config["ambient_pressure_psia"] * PSI_TO_PA, 1.0, mixture["specific_gas_constant_j_kg_k"], mixture["gamma"])
    sweep_design_flow = np.full_like(inflow, combined_design_flow)
    sweep_required_cda = sweep_design_flow / sweep_unit_capacity * M2_TO_MM2

    sensitivity_temperatures = list(np.linspace(thermal["sensitivity_min_temperature_k"], thermal["sensitivity_max_temperature_k"], 9))
    sensitivity_temperatures.append(relieving_temperature_k)
    if thermal["adiabatic_limiting_case"] is not None:
        sensitivity_temperatures.append(thermal["adiabatic_limiting_case"]["temperature_k"])
    sensitivity_rows = []
    for temperature_k in sorted(set(round(value, 8) for value in sensitivity_temperatures)):
        sensitivity_mixture = ideal_gas_mixture_properties(maximum_inflow, gox["mass_flow_kg_s"], temperature_k)
        sensitivity_capacity, sensitivity_choked = compressible_mass_flow(
            relief_pressure_psia * PSI_TO_PA, temperature_k,
            config["ambient_pressure_psia"] * PSI_TO_PA, 1.0,
            sensitivity_mixture["specific_gas_constant_j_kg_k"], sensitivity_mixture["gamma"],
        )
        sensitivity_cda = combined_design_flow / float(sensitivity_capacity) * M2_TO_MM2
        sensitivity_rows.append({
            "temperature_k": temperature_k,
            "required_cda_mm2": sensitivity_cda,
            "candidate_margin_mm2": config["orv"]["candidate_cda_mm2"] - sensitivity_cda,
            "relief_is_choked": bool(sensitivity_choked),
            "is_selected": math.isclose(temperature_k, relieving_temperature_k, abs_tol=1e-6),
            "is_adiabatic_limiting_case": thermal["adiabatic_limiting_case"] is not None and math.isclose(temperature_k, thermal["adiabatic_limiting_case"]["temperature_k"], abs_tol=1e-6),
        })
    return {
        **path, "gn2_properties": gn2, "relieving_gn2_properties": relieving_gn2, "mixture_properties": mixture, "gox": gox, "thermal": thermal,
        "pressure_psia": pressures_psia, "inflow_kg_s": inflow, "inflow_choked": inflow_choked,
        "pressure_matched_design_flow_kg_s": sweep_design_flow, "pressure_matched_required_cda_mm2": sweep_required_cda,
        "relief_choked_sweep": sweep_relief_choked,
        "maximum_gn2_inflow_kg_s": maximum_inflow, "analytical_maximum_gn2_inflow_kg_s": analytical_maximum_inflow,
        "analytical_inflow_difference_percent": analytical_difference_percent,
        "maximum_inflow_pressure_psia": maximum_inflow_pressure,
        "gn2_design_flow_kg_s": design_gn2, "combined_design_flow_kg_s": combined_design_flow,
        "governing_relieving_pressure_psia": relief_pressure_psia, "relief_is_choked": bool(relief_choked),
        "required_gn2_cda_mm2": required_gn2_cda, "required_cda_mm2": required_cda, "candidate_cda_mm2": candidate_cda,
        "required_diameter_mm": required_diameter,
        "candidate_capacity_kg_s": candidate_capacity, "cda_margin_mm2": candidate_cda - required_cda,
        "candidate_to_required_ratio": candidate_cda / required_cda,
        "preliminary_status": "PASS_ANALYTICAL_SCREENING_ONLY" if candidate_cda >= required_cda else "FAIL_ANALYTICAL_SCREENING",
        "temperature_sensitivity": sensitivity_rows,
        "temperature_sensitivity_cda_min_mm2": min(row["required_cda_mm2"] for row in sensitivity_rows),
        "temperature_sensitivity_cda_max_mm2": max(row["required_cda_mm2"] for row in sensitivity_rows),
        "bounding_method": f"simulation-derived maximum GN2 inflow combined with capacity at the configured {relief_pressure_psia:.3f} psia maximum relieving pressure",
    }


def _write_csv(path: Path, headers: list[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def write_outputs(config: dict[str, Any], result: dict[str, Any]) -> list[Path]:
    """Write the five-section summary plus detailed audit/sweep artifacts."""
    output = config["paths"]["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    # These belonged to the removed spreadsheet audit and structural
    # sensitivity load cases. Avoid leaving stale artifacts after a refactor run.
    for legacy_name in (
        "orv_relieving_pressure_scenarios.csv",
        "orv_spreadsheet_reference_audit.csv",
        "orv_load_case_sensitivity.png",
    ):
        (output / legacy_name).unlink(missing_ok=True)
    written: list[Path] = []
    summary = output / "orv_sizing_summary.csv"
    rows = [
        ("Requirements", "failure_scenario", config["failure"]["name"], "-", "Only modeled failure"),
        ("Requirements", "protected_system", config["tank"]["name"], "-", "LOX tank"),
        ("Requirements", "tank_meop", config["tank"]["meop_psia"], "psia", "Maximum expected operating pressure"),
        ("Requirements", "tank_mawp", config["tank"]["mawp_psia"], "psia", "Safety constraint; not operating target"),
        ("Requirements", "proposed_set_pressure", config["orv"]["proposed_set_pressure_psia"], "psia", "Not approved"),
        ("Requirements", "evaluated_relieving_pressure", result["governing_relieving_pressure_psia"], "psia", "Matches configured maximum relieving pressure"),
        ("Requirements", "maximum_relieving_pressure", config["orv"]["maximum_relieving_pressure_psia"], "psia", config["orv"]["provisional_marp_status"]),
        ("Requirements", "pressure_reference_convention", "absolute unless explicitly psig", "-", "Evaluated and maximum relieving pressures are both absolute"),
        ("Inputs", "gn2_source_pressure", max(config["source"]["pressure_psia"]), "psia", "hopper_config.json"),
        ("Inputs", "gn2_source_temperature", min(config["source"]["temperature_k"]), "K", "hopper_config.json"),
        ("Inputs", "failed_open_path_cda", result["path_cda_m2"] * M2_TO_MM2, "mm^2", "Derived SeriesCdA"),
        ("Inputs", "upstream_orifice_cda", config["failure"]["orifice_cda_mm2"], "mm^2", f"{config['failure']['failing_orifice_name']} from hopper_config.json"),
        ("Inputs", "obb_cda", config["failure"]["obb_cda_mm2"], "mm^2", "OBB from hopper_config.json"),
        ("Inputs", "simulation_maximum_gn2_inflow", result["maximum_gn2_inflow_kg_s"], "kg/s", config["failure"]["simulation_inflow_status"]),
        ("Inputs", "analytical_maximum_gn2_inflow", result["analytical_maximum_gn2_inflow_kg_s"], "kg/s", "Independent SeriesCdA screening estimate"),
        ("Inputs", "analytical_minus_simulation", result["analytical_inflow_difference_percent"], "%", "Does not replace governing simulation load"),
        ("Inputs", "orv_downstream_pressure", config["ambient_pressure_psia"], "psia", "Plumbing loss/backpressure neglected"),
        ("Inputs", "lox_temperature", config["tank"]["liquid_temperature_k"], "K", "hopper_config.json"),
        ("Inputs", "selected_relieving_temperature", result["thermal"]["selected_temperature_k"], "K", result["thermal"]["status"]),
        ("Inputs", "gox_generation", result["gox"]["mass_flow_kg_s"], "kg/s", result["gox"]["status"]),
        ("Inputs", "candidate_valve_cda", result["candidate_cda_mm2"], "mm^2", "Analytical effective area; not certified capacity"),
        ("Assumptions", "gn2_mass_flow_safety_factor", config["failure"]["gn2_mass_flow_safety_factor"], "-", "Applied once to GN2 only"),
        ("Assumptions", "mixture_model", "CoolProp component CP0 ideal-gas homogeneous mixture", "-", "Gas-only screening"),
        ("Assumptions", "analysis_model", "steady_state_gas_only", "-", "No transient or two-phase model"),
        ("Assumptions", "temperature_method", result["thermal"]["method"], "-", "Well-mixed energy-balance framework; selected result depends on thermal closure"),
        ("Assumptions", "missing_thermal_inputs", "; ".join(result["thermal"]["missing_inputs"]), "-", "See orv_thermal_balance.csv"),
        ("Assumptions", "bounding_method", result["bounding_method"], "-", "Explicit conservative bound"),
        ("Methodology", "maximum_gn2_inflow", "isentropic compressible flow through repository SeriesCdA", "-", "Choked/unchoked selected automatically"),
        ("Methodology", "combined_relief_load", f"{config['failure']['gn2_mass_flow_safety_factor']:.6g} * GN2 + GOX", "kg/s", "Safety factor is not applied twice"),
        ("Methodology", "required_cda", "combined design flow / mixture unit-area capacity", "-", "Uses inlet/outlet absolute pressures"),
        ("Methodology", "candidate_comparison", "candidate CdA - required CdA and ratio", "-", "Analytical screening only"),
        ("Results", "maximum_gn2_inflow", result["maximum_gn2_inflow_kg_s"], "kg/s", "Simulation-derived governing load"),
        ("Results", "combined_design_relief_flow", result["combined_design_flow_kg_s"], "kg/s", "Linked GN2 factor and GOX load"),
        ("Results", "governing_relieving_pressure", result["governing_relieving_pressure_psia"], "psia", "Matches provisional maximum; still pending approval"),
        ("Results", "relieving_temperature", result["thermal"]["selected_temperature_k"], "K", result["thermal"]["method"]),
        ("Results", "net_ullage_heat_transfer", result["thermal"]["net_heat_transfer_w"], "W", "Unresolved when no thermal closure is supplied"),
        ("Results", "implied_heat_transfer_at_selected_temperature", result["thermal"]["implied_heat_transfer_at_selected_temperature_w"], "W", "Heat into gas required to sustain selected temperature; negative means removal; not a predicted heat-transfer rate"),
        ("Results", "temperature_sensitivity_range", f"{result['thermal']['sensitivity_min_temperature_k']:.1f} to {result['thermal']['sensitivity_max_temperature_k']:.1f}", "K", "LOX-to-GN2 source-temperature bracket"),
        ("Results", "temperature_sensitivity_required_cda_range", f"{result['temperature_sensitivity_cda_min_mm2']:.4f} to {result['temperature_sensitivity_cda_max_mm2']:.4f}", "mm^2", "Detailed points in orv_temperature_sensitivity.csv"),
        ("Results", "gn2_only_required_orv_cda", result["required_gn2_cda_mm2"], "mm^2", "Regression comparison; linked safety factor on GN2"),
        ("Results", "required_orv_cda", result["required_cda_mm2"], "mm^2", "Preliminary analytical requirement"),
        ("Results", "required_equivalent_effective_diameter", result["required_diameter_mm"], "mm", "Derived from effective CdA"),
        ("Results", "candidate_orv_cda", result["candidate_cda_mm2"], "mm^2", "Not certified capacity"),
        ("Results", "candidate_theoretical_mass_flow_capacity", result["candidate_capacity_kg_s"], "kg/s", "Fully-open homogeneous-mixture screen"),
        ("Results", "candidate_cda_margin", result["cda_margin_mm2"], "mm^2", "Negative means shortfall"),
        ("Results", "candidate_to_required_ratio", result["candidate_to_required_ratio"], "-", "Ratio below one fails screening"),
        ("Results", "preliminary_sizing_status", result["preliminary_status"], "-", "Not approval or flight certification"),
        ("Results", "thermal_model_validation_status", "UNVALIDATED_PROVISIONAL_150_K", "-", "Thermal closure unavailable"),
        ("Results", "final_validation_status", "UNRESOLVED", "-", "Certified capacity, pressure approval, thermal closure, and two-phase applicability unresolved"),
    ]
    _write_csv(summary, ["section", "metric", "value", "unit", "note"], rows)
    written.append(summary)

    mapping = output / "orv_input_mapping.csv"
    headers = ["python_variable", "source_path", "units", "purpose"]
    _write_csv(mapping, headers, ([row[key] for key in headers] for row in config["input_mapping"]))
    written.append(mapping)
    unresolved = output / "orv_unresolved_inputs.csv"
    _write_csv(unresolved, ["input_requiring_engineering_confirmation"], ((item,) for item in config["unresolved_inputs"]))
    written.append(unresolved)
    conflicts = output / "orv_configuration_conflicts.csv"
    _write_csv(conflicts, ["configuration_conflict"], ((item,) for item in config["configuration_conflicts"]))
    written.append(conflicts)
    audit = output / "orv_input_audit.csv"
    audit_rows = [
        ("gn2_R", result["gn2_properties"]["specific_gas_constant_j_kg_k"], "J/(kg K)", "CoolProp calculated"),
        ("gn2_cp", result["gn2_properties"]["cp_j_kg_k"], "J/(kg K)", "CoolProp CP0"),
        ("gn2_cv", result["gn2_properties"]["cv_j_kg_k"], "J/(kg K)", "Derived cp-R"),
        ("gn2_gamma", result["gn2_properties"]["gamma"], "-", "Derived cp/cv"),
        ("relieving_gn2_gamma", result["relieving_gn2_properties"]["gamma"], "-", "CoolProp CP0 at relieving temperature"),
        ("mixture_R", result["mixture_properties"]["specific_gas_constant_j_kg_k"], "J/(kg K)", "Mass-weighted ideal mixture"),
        ("mixture_gamma", result["mixture_properties"]["gamma"], "-", "Derived cp/cv"),
        ("relief_is_choked", int(result["relief_is_choked"]), "0_or_1", "Calculated"),
        ("candidate_capacity", result["candidate_capacity_kg_s"], "kg/s", "Analytical mixture capacity"),
    ]
    _write_csv(audit, ["derived_parameter", "value", "unit", "note"], audit_rows)
    written.append(audit)

    thermal_output = output / "orv_thermal_balance.csv"
    adiabatic = result["thermal"]["adiabatic_limiting_case"]
    calculated = result["thermal"]["calculated"]
    thermal_rows = [
        ("temperature_selection_method", result["thermal"]["method"], "-", result["thermal"]["status"]),
        ("selected_relieving_temperature", result["thermal"]["selected_temperature_k"], "K", "Used for ORV CdA"),
        ("gn2_physical_mass_flow", result["maximum_gn2_inflow_kg_s"], "kg/s", "Safety factor intentionally excluded from energy balance"),
        ("gn2_inlet_temperature", min(config["source"]["temperature_k"]), "K", "COPV stagnation temperature; throttling treated as isenthalpic"),
        ("gox_physical_mass_flow", result["gox"]["mass_flow_kg_s"], "kg/s", result["gox"]["status"]),
        ("gox_inlet_temperature", result["thermal"]["gox_inlet_temperature_k"], "K", "Saturated-vapor state; latent heat not counted again"),
        ("lox_temperature", config["tank"]["liquid_temperature_k"], "K", "hopper_config.json"),
        ("liquid_to_ullage_ua", config["thermal"]["liquid_to_ullage_ua_w_k"], "W/K", "Missing" if config["thermal"]["liquid_to_ullage_ua_w_k"] is None else "Input"),
        ("wall_to_ullage_ua", config["thermal"]["wall_to_ullage_ua_w_k"], "W/K", "Missing" if config["thermal"]["wall_to_ullage_ua_w_k"] is None else "Input"),
        ("tank_wall_temperature", config["thermal"]["tank_wall_temperature_k"], "K", "Missing" if config["thermal"]["tank_wall_temperature_k"] is None else "Input"),
        ("net_external_sensible_heat_rate", config["thermal"]["net_external_sensible_heat_rate_w"], "W", "Missing; must exclude vaporization energy already represented by GOX flow"),
        ("adiabatic_well_mixed_limiting_temperature", None if adiabatic is None else adiabatic["temperature_k"], "K", "Calculated limiting case; not selected without thermal justification"),
        ("calculated_relieving_temperature", None if calculated is None else calculated["temperature_k"], "K", "Blank when thermal closure is unavailable"),
        ("calculated_net_heat_transfer", None if calculated is None else calculated["net_heat_transfer_w"], "W", "Positive into ullage gas"),
        ("implied_heat_transfer_at_selected_temperature", result["thermal"]["implied_heat_transfer_at_selected_temperature_w"], "W", "Required by selected temperature, not supported by a UA/wall model"),
        ("energy_balance_residual", None if calculated is None else calculated["energy_residual_w"], "W", "Numerical convergence check"),
        ("missing_inputs", "; ".join(result["thermal"]["missing_inputs"]), "-", "Blocks calculated temperature selection"),
    ]
    _write_csv(thermal_output, ["thermal_metric", "value", "unit", "note"], thermal_rows)
    written.append(thermal_output)

    sensitivity_output = output / "orv_temperature_sensitivity.csv"
    _write_csv(
        sensitivity_output,
        ["relieving_temperature_k", "required_orv_cda_mm2", "candidate_cda_margin_mm2", "relief_is_choked", "is_selected_temperature", "is_adiabatic_limiting_case"],
        ((row["temperature_k"], row["required_cda_mm2"], row["candidate_margin_mm2"], int(row["relief_is_choked"]), int(row["is_selected"]), int(row["is_adiabatic_limiting_case"])) for row in result["temperature_sensitivity"]),
    )
    written.append(sensitivity_output)
    sweep = output / "orv_pressure_sweep.csv"
    _write_csv(sweep, ["tank_pressure_psia", "gn2_inflow_kg_s", "inflow_is_choked", "pressure_matched_design_flow_kg_s", "pressure_matched_required_cda_mm2", "relief_is_choked"], zip(result["pressure_psia"], result["inflow_kg_s"], result["inflow_choked"].astype(int), result["pressure_matched_design_flow_kg_s"], result["pressure_matched_required_cda_mm2"], result["relief_choked_sweep"].astype(int)))
    written.append(sweep)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(result["pressure_psia"], result["inflow_kg_s"], label="Analytical failed-open GN2 inflow")
    axis.axhline(result["maximum_gn2_inflow_kg_s"], color="tab:red", label="Simulation maximum GN2 (governing)")
    axis.plot(result["pressure_psia"], result["pressure_matched_design_flow_kg_s"], "--", label=f"Governing {config['failure']['gn2_mass_flow_safety_factor']:.3g} × GN2 + GOX")
    axis.axvline(result["governing_relieving_pressure_psia"], color="black", linestyle=":", label="Exploratory relief pressure")
    axis.set(xlabel="Tank pressure [psia]", ylabel="Mass flow [kg/s]", title="Oxidizer B/B fail-open ORV flow screening")
    axis.grid(True); axis.legend()
    flow_plot = output / "orv_flow_comparison.png"
    fig.savefig(flow_plot, dpi=200); plt.close(fig); written.append(flow_plot)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(result["pressure_psia"], result["pressure_matched_required_cda_mm2"], label="Pressure-matched required CdA")
    axis.axhline(result["required_cda_mm2"], color="black", linestyle="--", label="Conservative bounding required CdA")
    axis.axhline(result["candidate_cda_mm2"], color="tab:orange", linestyle=":", label="Candidate CdA")
    axis.set(xlabel="Tank/relieving pressure [psia]", ylabel="Effective CdA [mm²]", title="Preliminary ORV effective area — not certified")
    axis.grid(True); axis.legend()
    cda_plot = output / "orv_required_cda.png"
    fig.savefig(cda_plot, dpi=200); plt.close(fig); written.append(cda_plot)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    temperatures = [row["temperature_k"] for row in result["temperature_sensitivity"]]
    required_areas = [row["required_cda_mm2"] for row in result["temperature_sensitivity"]]
    axis.plot(temperatures, required_areas, marker="o", label="Combined GN2/GOX required CdA")
    axis.axhline(result["candidate_cda_mm2"], color="tab:orange", linestyle="--", label="Candidate CdA")
    axis.axvline(result["thermal"]["selected_temperature_k"], color="black", linestyle=":", label="Selected provisional temperature")
    axis.set(xlabel="Relieving temperature [K]", ylabel="Required effective CdA [mm²]", title="ORV temperature sensitivity — gas-only screening")
    axis.grid(True); axis.legend()
    temperature_plot = output / "orv_required_cda_vs_temperature.png"
    fig.savefig(temperature_plot, dpi=200); plt.close(fig); written.append(temperature_plot)
    return written


def print_report(config: dict[str, Any], result: dict[str, Any], written: list[Path]) -> None:
    print("\nPRELIMINARY OXIDIZER B/B FAIL-OPEN ORV SIZING — NOT CERTIFIED")
    print(f"Maximum GN2 inflow: {result['maximum_gn2_inflow_kg_s']:.7f} kg/s")
    print(f"Provisional GOX generation: {result['gox']['mass_flow_kg_s']:.6f} kg/s ({result['gox']['status']})")
    print(f"Combined design flow ({config['failure']['gn2_mass_flow_safety_factor']:.6g} * GN2 + GOX): {result['combined_design_flow_kg_s']:.7f} kg/s")
    print(f"Evaluated/maximum relieving pressure: {result['governing_relieving_pressure_psia']:.2f} psia")
    print(f"Relieving temperature: {result['thermal']['selected_temperature_k']:.2f} K ({result['thermal']['method']})")
    if result["thermal"]["adiabatic_limiting_case"] is not None:
        print(f"Adiabatic well-mixed limiting case: {result['thermal']['adiabatic_limiting_case']['temperature_k']:.2f} K (not selected)")
    print(f"Temperature sensitivity CdA: {result['temperature_sensitivity_cda_min_mm2']:.4f}–{result['temperature_sensitivity_cda_max_mm2']:.4f} mm^2")
    print(f"GN2-only required CdA: {result['required_gn2_cda_mm2']:.4f} mm^2")
    print(f"Required / candidate CdA: {result['required_cda_mm2']:.4f} / {result['candidate_cda_mm2']:.4f} mm^2")
    print(f"Margin / ratio: {result['cda_margin_mm2']:.4f} mm^2 / {result['candidate_to_required_ratio']:.4f}")
    print(f"Status: {result['preliminary_status']} (analytical only)")
    print(f"Bounding method: {result['bounding_method']}")
    print("\nEngineering confirmation required:")
    for item in config["unresolved_inputs"]:
        print(f"  - {item}")
    print("\nWritten files:")
    for path in written:
        print(f"  - {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_HOPPER_CONFIG, help="Read-only Hopper JSON path")
    arguments = parser.parse_args()
    config = load_orv_config(arguments.config)
    warnings.warn("SeriesCdA and ideal-gas homogeneous relief are preliminary screening approximations", stacklevel=2)
    result = calculate_sizing(config)
    written = write_outputs(config, result)
    print_report(config, result, written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Preliminary steady-state analytical sizing of the Hopper fuel relief valve.

Equations
---------
The non-inertial gas line is represented exactly as documented by the Hopper
model::

    f = [-2 log10(epsilon / (3.7 D))]^-2
    CdA_line = A / sqrt(1 + f L / D)

Configured restrictions are collapsed with the FullFlow SeriesCdA convention::

    CdA_eq = [sum_i(1 / CdA_i^2)]^-1/2

Ideal-gas isentropic orifice flow is evaluated with absolute pressures.  For
downstream/upstream pressure ratio r at or below the critical ratio, flow is::

    mdot = CdA P0 sqrt(gamma / (R T0))
           [2 / (gamma + 1)]^[(gamma + 1) / (2(gamma - 1))]

Above the critical ratio, the corresponding unchoked isentropic expression is
used.  The theoretical equality threshold and margin-adjusted design
requirement are::

    mdot_FRV_theoretical = max(mdot_FBB_fail_open)
    mdot_FRV_design = (1 + M) max(mdot_FBB_fail_open),  M > 0

These equations are screening tools.  They are not a manufacturer-certified
capacity method and do not establish final hardware suitability.
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

from sims.hopper_fullflow.isolated_sims.FRV_sizing.frv_config import DEFAULT_HOPPER_CONFIG, load_frv_config


PSI_TO_PA = 6894.757293168
INCH_TO_M = 0.0254
MM2_TO_M2 = 1e-6
M2_TO_MM2 = 1e6
M_TO_MM = 1e3


def series_equivalent_cda(areas_m2: Iterable[float]) -> float:
    """Return the repository SeriesCdA effective area for restrictions in series."""
    areas = np.asarray(tuple(areas_m2), dtype=float)
    if areas.ndim != 1 or areas.size == 0 or np.any(~np.isfinite(areas)) or np.any(areas <= 0):
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
    """Ideal-gas isentropic orifice flow and choking state.

    Every pressure is absolute.  ``effective_cda_m2`` already includes the
    discharge coefficient.  Inputs are explicitly broadcast to prevent silent
    row/column behavior differences between MATLAB and NumPy.
    """
    upstream, temperature, downstream, cda = np.broadcast_arrays(
        np.asarray(upstream_pressure_pa, dtype=float),
        np.asarray(upstream_temperature_k, dtype=float),
        np.asarray(downstream_pressure_pa, dtype=float),
        np.asarray(effective_cda_m2, dtype=float),
    )
    if np.any(~np.isfinite(upstream)) or np.any(upstream <= 0):
        raise ValueError("Upstream absolute pressure must be finite and positive")
    if np.any(~np.isfinite(downstream)) or np.any(downstream <= 0):
        raise ValueError("Downstream absolute pressure must be finite and positive")
    if np.any(~np.isfinite(temperature)) or np.any(temperature <= 0):
        raise ValueError("Upstream absolute temperature must be finite and positive")
    if np.any(~np.isfinite(cda)) or np.any(cda <= 0):
        raise ValueError("Effective CdA must be finite and positive")

    ratio = downstream / upstream
    critical_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    forward = upstream > downstream
    choked = forward & (ratio <= critical_ratio)
    unchoked = forward & ~choked
    mass_flow = np.zeros(upstream.shape, dtype=float)

    choked_factor = np.sqrt(gamma / (gas_constant_j_kg_k * temperature)) * (
        2.0 / (gamma + 1.0)
    ) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    mass_flow[choked] = cda[choked] * upstream[choked] * choked_factor[choked]

    unchoked_term = ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma)
    unchoked_factor = np.sqrt(
        np.maximum(
            0.0,
            2.0 * gamma / (gas_constant_j_kg_k * temperature * (gamma - 1.0))
            * unchoked_term,
        )
    )
    mass_flow[unchoked] = cda[unchoked] * upstream[unchoked] * unchoked_factor[unchoked]
    return mass_flow, choked


def linear_lift_fraction(pressure_psia: np.ndarray, set_pressure_psia: float, full_lift_pressure_psia: float) -> np.ndarray:
    """Exploratory, non-certified linear opening fraction."""
    if full_lift_pressure_psia <= set_pressure_psia:
        raise ValueError("Exploratory full-lift pressure must exceed proposed set pressure")
    return np.clip(
        (pressure_psia - set_pressure_psia) / (full_lift_pressure_psia - set_pressure_psia),
        0.0,
        1.0,
    )


def calculate(config: dict[str, Any]) -> dict[str, Any]:
    """Calculate the preliminary FRV pressure sweep without writing files."""
    source = config["source"]
    tank = config["tank"]
    failure = config["failure"]
    frv = config["frv"]
    pressurant = config["pressurant"]

    ambient_pa = config["ambient_pressure_psia"] * PSI_TO_PA
    source_pressure_pa = max(source["pressure_psia"]) * PSI_TO_PA
    source_temperature_k = min(source["temperature_k"])
    ullage_temperature_k = tank["ullage_temperature_k"]
    gas_constant = pressurant["specific_gas_constant_j_kg_k"]
    gamma = pressurant["specific_heat_ratio"]

    line = failure["pressurant_line"]
    diameter_in = line.get("inner_diameter_in")
    if diameter_in is None:
        diameter_in = line["outer_diameter_in"] - 2.0 * line["wall_thickness_in"]
    diameter_m = diameter_in * INCH_TO_M
    length_m = line["length_in"] * INCH_TO_M
    roughness_m = line["roughness_in"] * INCH_TO_M
    area_m2 = math.pi * diameter_m**2 / 4.0
    friction_factor = (-2.0 * math.log10(roughness_m / (3.7 * diameter_m))) ** -2.0
    line_cda_m2 = area_m2 / math.sqrt(1.0 + friction_factor * length_m / diameter_m)

    bottle_outlet_cda_m2 = sum(source["outlet_cda_mm2"]) * MM2_TO_M2
    path_cda_m2 = series_equivalent_cda(
        (
            bottle_outlet_cda_m2,
            line_cda_m2,
            failure["fbb_cda_mm2"] * MM2_TO_M2,
            failure["orifice_cda_mm2"] * MM2_TO_M2,
            failure["check_valve_cda_mm2"] * MM2_TO_M2,
        )
    )

    set_pressure_psia = frv["proposed_set_pressure_psia"]
    full_lift_pressure_psia = set_pressure_psia + frv["full_lift_overpressure_fraction"] * (
        set_pressure_psia - config["ambient_pressure_psia"]
    )
    maximum_relief_known = tank["maximum_allowable_relieving_pressure_psia"] is not None
    if maximum_relief_known:
        evaluation_pressure_psia = tank["maximum_allowable_relieving_pressure_psia"]
        evaluation_basis = "maximum_allowable_relieving_pressure"
    else:
        evaluation_pressure_psia = full_lift_pressure_psia
        evaluation_basis = "exploratory_full_lift_pressure_MAWP_and_max_relief_unknown"
    if evaluation_pressure_psia < full_lift_pressure_psia:
        raise ValueError(
            "The relieving-pressure evaluation point is below the assumed full-open pressure"
        )
    if source_pressure_pa <= evaluation_pressure_psia * PSI_TO_PA:
        raise ValueError("COPV source pressure must exceed the analytical evaluation pressure")

    pressure_psia = np.linspace(
        tank["meop_psia"], evaluation_pressure_psia, config["analysis"]["pressure_sweep_points"]
    )
    breakpoints = [set_pressure_psia, evaluation_pressure_psia]
    if full_lift_pressure_psia <= evaluation_pressure_psia:
        breakpoints.append(full_lift_pressure_psia)
    pressure_psia = np.unique(np.concatenate((pressure_psia, np.asarray(breakpoints))))
    pressure_pa = pressure_psia * PSI_TO_PA

    inflow_kg_s, inflow_choked = compressible_mass_flow(
        source_pressure_pa, source_temperature_k, pressure_pa, path_cda_m2, gas_constant, gamma
    )
    maximum_index = int(np.argmax(inflow_kg_s))
    analytical_maximum_inflow_kg_s = float(inflow_kg_s[maximum_index])
    governing_inflow_kg_s = failure["simulation_maximum_gn2_inflow_kg_s"]
    inflow_difference_percent = 100.0 * (
        analytical_maximum_inflow_kg_s - governing_inflow_kg_s
    ) / governing_inflow_kg_s
    sizing_margin_fraction = failure["sizing_margin_fraction"]
    theoretical_flow_kg_s = governing_inflow_kg_s
    design_flow_kg_s = (1.0 + sizing_margin_fraction) * governing_inflow_kg_s
    local_design_flow_kg_s = np.full_like(inflow_kg_s, design_flow_kg_s)

    unit_area_capacity_kg_s_m2, outflow_choked = compressible_mass_flow(
        pressure_pa, ullage_temperature_k, ambient_pa, 1.0, gas_constant, gamma
    )
    theoretical_cda_m2 = theoretical_flow_kg_s / unit_area_capacity_kg_s_m2
    design_cda_m2 = design_flow_kg_s / unit_area_capacity_kg_s_m2
    theoretical_cda_mm2 = theoretical_cda_m2 * M2_TO_MM2
    design_cda_mm2 = design_cda_m2 * M2_TO_MM2
    design_diameter_mm = np.sqrt(4.0 * design_cda_m2 / math.pi) * M_TO_MM

    configured_cda_m2 = frv["configured_full_open_cda_mm2"] * MM2_TO_M2
    configured_envelope_kg_s, _ = compressible_mass_flow(
        pressure_pa, ullage_temperature_k, ambient_pa, configured_cda_m2, gas_constant, gamma
    )
    lift_fraction = linear_lift_fraction(pressure_psia, set_pressure_psia, full_lift_pressure_psia)
    exploratory_outflow_kg_s = lift_fraction * configured_envelope_kg_s
    candidate_to_incoming_ratio = exploratory_outflow_kg_s / governing_inflow_kg_s
    candidate_to_design_ratio = exploratory_outflow_kg_s / design_flow_kg_s
    exploratory_margin_percent = 100.0 * (candidate_to_design_ratio - 1.0)
    evaluation_index = int(np.argmin(np.abs(pressure_psia - evaluation_pressure_psia)))

    return {
        "pressure_psia": pressure_psia,
        "inflow_kg_s": inflow_kg_s,
        "inflow_choked": inflow_choked,
        "local_design_flow_kg_s": local_design_flow_kg_s,
        "theoretical_cda_mm2": theoretical_cda_mm2,
        "design_cda_mm2": design_cda_mm2,
        "design_diameter_mm": design_diameter_mm,
        "configured_envelope_kg_s": configured_envelope_kg_s,
        "outflow_choked": outflow_choked,
        "lift_fraction": lift_fraction,
        "exploratory_outflow_kg_s": exploratory_outflow_kg_s,
        "candidate_to_incoming_ratio": candidate_to_incoming_ratio,
        "candidate_to_design_ratio": candidate_to_design_ratio,
        "exploratory_margin_percent": exploratory_margin_percent,
        "line_cda_mm2": line_cda_m2 * M2_TO_MM2,
        "path_cda_mm2": path_cda_m2 * M2_TO_MM2,
        "full_lift_pressure_psia": full_lift_pressure_psia,
        "evaluation_pressure_psia": evaluation_pressure_psia,
        "evaluation_basis": evaluation_basis,
        "maximum_inflow_kg_s": governing_inflow_kg_s,
        "analytical_maximum_inflow_kg_s": analytical_maximum_inflow_kg_s,
        "analytical_inflow_difference_percent": inflow_difference_percent,
        "maximum_inflow_pressure_psia": float(pressure_psia[maximum_index]),
        "theoretical_flow_kg_s": theoretical_flow_kg_s,
        "design_flow_kg_s": design_flow_kg_s,
        "theoretical_cda_at_evaluation_mm2": float(theoretical_cda_mm2[evaluation_index]),
        "design_cda_at_evaluation_mm2": float(design_cda_mm2[evaluation_index]),
        "design_diameter_at_evaluation_mm": float(design_diameter_mm[evaluation_index]),
        "candidate_capacity_at_evaluation_kg_s": float(exploratory_outflow_kg_s[evaluation_index]),
        "candidate_to_incoming_ratio_at_evaluation": float(candidate_to_incoming_ratio[evaluation_index]),
        "candidate_to_design_ratio_at_evaluation": float(candidate_to_design_ratio[evaluation_index]),
        "exploratory_margin_at_evaluation_percent": float(exploratory_margin_percent[evaluation_index]),
        "candidate_satisfies_design_requirement": bool(exploratory_outflow_kg_s[evaluation_index] >= design_flow_kg_s),
        "configured_hardware_crack_equivalent_psia": frv["configured_hardware_crack_pressure_psig"] + config["ambient_pressure_psia"],
        "maximum_relief_known": maximum_relief_known,
    }


def _write_csv(path: Path, headers: list[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def write_outputs(config: dict[str, Any], result: dict[str, Any]) -> list[Path]:
    """Write traceable CSVs and non-interactive plots to Outputs/FRV_outputs."""
    output = config["paths"]["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    sweep_path = output / "frv_pressure_sweep.csv"
    _write_csv(
        sweep_path,
        [
            "tank_pressure_psia", "fbb_fail_open_inflow_kg_s", "inflow_is_choked",
            "design_flow_at_local_inflow_kg_s", "theoretical_effective_cda_mm2",
            "design_effective_cda_mm2", "design_equivalent_diameter_mm",
            "configured_full_open_analytical_envelope_kg_s",
            "frv_outflow_is_choked", "exploratory_lift_fraction",
            "candidate_frv_outflow_kg_s", "candidate_to_incoming_flow_ratio",
            "candidate_to_design_flow_ratio", "candidate_margin_percent",
        ],
        zip(
            result["pressure_psia"], result["inflow_kg_s"], result["inflow_choked"].astype(int),
            result["local_design_flow_kg_s"], result["theoretical_cda_mm2"], result["design_cda_mm2"],
            result["design_diameter_mm"],
            result["configured_envelope_kg_s"], result["outflow_choked"].astype(int), result["lift_fraction"],
            result["exploratory_outflow_kg_s"], result["candidate_to_incoming_ratio"],
            result["candidate_to_design_ratio"], result["exploratory_margin_percent"],
        ),
    )
    written.append(sweep_path)

    summary_rows = [
        ("Requirements", "failure_scenario", config["failure"]["name"], "-", "FBB fails open"),
        ("Requirements", "protected_system", config["tank"]["name"], "-", "Fuel tank"),
        ("Requirements", "tank_meop", config["tank"]["meop_psia"], "psia", "hopper_config.json"),
        ("Requirements", "tank_mawp", config["tank"]["mawp_psia"], "psia", config["tank"]["mawp_status"]),
        ("Requirements", "set_pressure", config["frv"]["proposed_set_pressure_psia"], "psia", config["frv"]["set_pressure_status"]),
        ("Requirements", "evaluated_relieving_pressure", result["evaluation_pressure_psia"], "psia", "Configured relief full-open pressure; not approved MARP"),
        ("Requirements", "pressure_reference_convention", "absolute unless explicitly psig", "-", "Set/full-open gauge pressures are converted with configured ambient pressure"),
        ("Inputs", "simulation_maximum_gn2_inflow", result["maximum_inflow_kg_s"], "kg/s", config["failure"]["simulation_inflow_status"]),
        ("Inputs", "analytical_maximum_gn2_inflow", result["analytical_maximum_inflow_kg_s"], "kg/s", "Independent SeriesCdA screening estimate"),
        ("Inputs", "analytical_minus_simulation", result["analytical_inflow_difference_percent"], "%", "Difference is not used to overwrite governing load"),
        ("Inputs", "analytical_inflow_discrepancy_flag", "SIGNIFICANT" if abs(result["analytical_inflow_difference_percent"]) > 10.0 else "WITHIN_10_PERCENT", "-", "Likely reflects SeriesCdA/multistage-model limitations or differing simulation conditions"),
        ("Inputs", "upstream_orifice_cda", config["failure"]["orifice_cda_mm2"], "mm^2", f"{config['failure']['orifice_name']} from hopper_config.json"),
        ("Inputs", "fbb_cda", config["failure"]["fbb_cda_mm2"], "mm^2", "FBB from hopper_config.json"),
        ("Inputs", "combined_path_effective_cda", result["path_cda_mm2"], "mm^2", f"Includes outlets, line, FBB, {config['failure']['orifice_name']}, PRFCK; SeriesCdA approximation"),
        ("Inputs", "gn2_relieving_temperature", config["tank"]["ullage_temperature_k"], "K", "Fuel ullage assumption"),
        ("Inputs", "atmospheric_discharge_pressure", config["ambient_pressure_psia"], "psia", "Discharge losses neglected"),
        ("Inputs", "candidate_valve_cda", config["frv"]["configured_full_open_cda_mm2"], "mm^2", "Configured theoretical fully-open CdA"),
        ("Assumptions", "mass_flow_factor", 1.0 + config["failure"]["sizing_margin_fraction"], "-", "Applied exactly once"),
        ("Assumptions", "gas_property_model", config["pressurant"]["model"], "-", "Constant-gamma GN2 screening"),
        ("Assumptions", "steady_state_limitations", "SeriesCdA and fully-open capacity are screening approximations", "-", "Not certified multistage or valve-lift performance"),
        ("Methodology", "governing_relief_load", f"{1.0 + config['failure']['sizing_margin_fraction']:.6g} * simulation maximum GN2 inflow", "-", "Analytical inflow retained only as validation"),
        ("Methodology", "compressible_flow", "isentropic choked/unchoked", "-", "Absolute pressures"),
        ("Results", "maximum_governing_inflow", result["maximum_inflow_kg_s"], "kg/s", config["failure"]["simulation_inflow_status"]),
        ("Results", "required_design_relief_flow", result["design_flow_kg_s"], "kg/s", "Linked maximum inflow times linked safety factor"),
        ("Results", "required_effective_cda", result["design_cda_at_evaluation_mm2"], "mm^2", "At evaluated relieving pressure"),
        ("Results", "equivalent_effective_diameter", result["design_diameter_at_evaluation_mm"], "mm", "From required CdA"),
        ("Results", "candidate_effective_cda", config["frv"]["configured_full_open_cda_mm2"], "mm^2", "Not certified"),
        ("Results", "candidate_theoretical_mass_flow_capacity", result["candidate_capacity_at_evaluation_kg_s"], "kg/s", "Exploratory full-lift envelope"),
        ("Results", "capacity_margin", result["exploratory_margin_at_evaluation_percent"], "%", "Against design flow"),
        ("Results", "candidate_screening", "PASS" if result["candidate_satisfies_design_requirement"] else "FAIL", "-", "Analytical only"),
        ("Results", "final_validation_status", "UNRESOLVED", "-", "Certified capacity, approved relieving pressure, and inlet/outlet losses unavailable"),
    ]
    summary_path = output / "frv_sizing_summary.csv"
    _write_csv(summary_path, ["section", "metric", "value", "unit", "note"], summary_rows)
    written.append(summary_path)

    missing_path = output / "frv_unresolved_inputs.csv"
    _write_csv(missing_path, ["unresolved_engineering_input"], ((item,) for item in config["missing_inputs"]))
    written.append(missing_path)

    mapping_path = output / "frv_input_mapping.csv"
    mapping_headers = ["python_variable", "json_source_path", "source_units", "purpose"]
    _write_csv(mapping_path, mapping_headers, ([row[key] for key in mapping_headers] for row in config["input_mapping"]))
    written.append(mapping_path)

    pressure = result["pressure_psia"]
    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(pressure, result["inflow_kg_s"], linewidth=2, label="FBB fail-open analytical inflow")
    axis.axhline(result["maximum_inflow_kg_s"], color="tab:red", linewidth=2, label="Simulation maximum (governing)")
    axis.plot(pressure, result["local_design_flow_kg_s"], "--", linewidth=2, label="Design flow with margin")
    axis.plot(pressure, result["configured_envelope_kg_s"], linewidth=2, label="Configured-CdA fully-open envelope")
    axis.plot(pressure, result["exploratory_outflow_kg_s"], ":", linewidth=2, label="Exploratory linear-lift outflow")
    axis.axvline(config["frv"]["proposed_set_pressure_psia"], color="black", linestyle=":", label="Proposed set")
    axis.axvline(result["full_lift_pressure_psia"], color="gray", linestyle=":", label="Exploratory full lift")
    axis.set(xlabel="Fuel-tank absolute pressure [psia]", ylabel="GN2 mass flow [kg/s]", title="Preliminary analytical FRV flow comparison — not certified capacity")
    axis.grid(True)
    axis.legend(loc="best")
    flow_plot = output / "frv_flow_comparison.png"
    fig.savefig(flow_plot, dpi=200)
    plt.close(fig)
    written.append(flow_plot)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(pressure, result["theoretical_cda_mm2"], linewidth=2, label="Minimum theoretical CdA")
    axis.plot(pressure, result["design_cda_mm2"], "--", linewidth=2, label="Design CdA with margin")
    axis.axhline(config["frv"]["configured_full_open_cda_mm2"], color="tab:orange", linestyle="--", label="Configured analytical CdA")
    axis.axvline(config["frv"]["proposed_set_pressure_psia"], color="black", linestyle=":", label="Proposed set")
    axis.axvline(result["full_lift_pressure_psia"], color="gray", linestyle=":", label="Exploratory full lift")
    axis.set(xlabel="Fuel-tank absolute pressure [psia]", ylabel="Required analytical effective CdA [mm²]", title="Preliminary required effective area — final sizing not established")
    axis.grid(True)
    axis.legend(loc="best")
    cda_plot = output / "frv_required_cda.png"
    fig.savefig(cda_plot, dpi=200)
    plt.close(fig)
    written.append(cda_plot)
    return written


def print_report(config: dict[str, Any], result: dict[str, Any], written: list[Path]) -> None:
    mawp_text = (
        "UNVERIFIED / NOT USED"
        if config["tank"]["mawp_psia"] is None
        else f"{config['tank']['mawp_psia']:.1f} psia"
    )
    relief_text = (
        "UNVERIFIED / NOT USED"
        if not result["maximum_relief_known"]
        else f"{config['tank']['maximum_allowable_relieving_pressure_psia']:.1f} psia"
    )
    print("\nPRELIMINARY FRV ANALYTICAL SIZING — NOT CERTIFIED CAPACITY")
    print(f"Hopper input source: {config['paths']['hopper_config']}")
    print("FRV sizing inputs source: hopper_config.json:relief_sizing.frv")
    print(f"Failure scenario: {config['failure']['name']}")
    print(f"Fuel-tank MEOP: {config['tank']['meop_psia']:.1f} psia ({config['source_paths']['tank_meop_psia']})")
    print(f"Proposed FRV set pressure: {config['frv']['proposed_set_pressure_psia']:.1f} psia ({config['frv']['set_pressure_status']})")
    print(f"Existing hardware crack input: {config['frv']['configured_hardware_crack_pressure_psig']:.1f} psig = {result['configured_hardware_crack_equivalent_psia']:.1f} psia at configured ambient")
    print(f"Fuel-tank MAWP: {mawp_text} ({config['tank']['mawp_status']})")
    print(f"Assumed permitted accumulation: {config['tank']['permitted_pressure_accumulation_psi']:.1f} psi above set")
    print(f"Maximum allowable relieving pressure: {relief_text} ({config['tank']['relieving_pressure_status']})")
    print(f"Analytical evaluation pressure: {result['evaluation_pressure_psia']:.1f} psia ({result['evaluation_basis']})")
    print(f"FBB failure path: COPV outlets -> Fuel Pressurant Line -> FBB -> {config['failure']['orifice_name']} -> PRFCK -> Fuel Ullage")
    print(f"Regulator treatment: {config['context']['regulator_name']} is a purge regulator and is excluded from the FBB path")
    print(f"Fuel pressurant equivalent CdA: {result['path_cda_mm2']:.4f} mm^2")
    print(f"Maximum FBB fail-open inflow: {result['maximum_inflow_kg_s']:.6f} kg/s")
    print(f"Positive sizing margin M: {config['failure']['sizing_margin_fraction']:.3f} ({100.0 * config['failure']['sizing_margin_fraction']:.1f} %)")
    print(f"FRV design outgoing flow: {result['design_flow_kg_s']:.6f} kg/s")
    print(f"Minimum theoretical FRV CdA (no margin): {result['theoretical_cda_at_evaluation_mm2']:.4f} mm^2")
    print(f"Design FRV CdA including margin: {result['design_cda_at_evaluation_mm2']:.4f} mm^2")
    print(f"Design equivalent effective diameter: {result['design_diameter_at_evaluation_mm']:.4f} mm")
    print(f"Candidate valve CdA: {config['frv']['configured_full_open_cda_mm2']:.4f} mm^2")
    print(f"Candidate discharge capacity: {result['candidate_capacity_at_evaluation_kg_s']:.6f} kg/s")
    print(f"Candidate outgoing/incoming flow ratio: {result['candidate_to_incoming_ratio_at_evaluation']:.4f}")
    print(f"Candidate satisfies steady-state design requirement: {'YES' if result['candidate_satisfies_design_requirement'] else 'NO'}")
    print("FINAL SIZING CONCLUSION: NOT ESTABLISHED (MARP/PPA and certified capacity unapproved).")
    print("MANUFACTURER-CERTIFIED VALVE CAPACITY: NOT AVAILABLE.")
    print("MODEL LIMIT: constant-gamma ideal gas and SeriesCdA are preliminary approximations.")
    print("SIMULTANEOUS FBB/OBB: requires a shared-source transient evaluation.")
    print("\nUnresolved engineering inputs:")
    for item in config["missing_inputs"]:
        print(f"  - {item}")
    print("\nWritten files:")
    for path in written:
        print(f"  - {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_HOPPER_CONFIG, help="Hopper JSON input (read-only)")
    arguments = parser.parse_args()
    config = load_frv_config(arguments.config)
    if config["source"]["model"] != "lumped":
        warnings.warn("Non-lumped COPVs require a shared-header solution; max-P/min-T is screening only", stacklevel=2)
    if config["tank"]["mawp_psia"] is None:
        warnings.warn("Fuel-tank MAWP value/units are unverified; final sizing is blocked", stacklevel=2)
    if config["tank"]["maximum_allowable_relieving_pressure_psia"] is None:
        warnings.warn("Maximum allowable relieving pressure is unverified; final sizing is blocked", stacklevel=2)
    warnings.warn("Constant-gamma ideal-gas flow and SeriesCdA are preliminary screening approximations", stacklevel=2)
    warnings.warn("The linear lift curve is exploratory and not manufacturer-certified", stacklevel=2)
    result = calculate(config)
    written = write_outputs(config, result)
    print_report(config, result, written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

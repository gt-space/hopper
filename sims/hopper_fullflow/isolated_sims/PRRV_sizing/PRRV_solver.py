"""High-pressure GN2 PRRV known-flow and unresolved-flow trade screening."""

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

try:
    from .PRRV_config import DEFAULT_HOPPER_CONFIG, load_prrv_inputs
except ImportError:  # Direct execution of a copied two-file handoff.
    from PRRV_config import DEFAULT_HOPPER_CONFIG, load_prrv_inputs

PSI_TO_PA = 6894.757293168
MM2_TO_M2 = 1e-6
M2_TO_MM2 = 1e6


def pressure_to_psia(value_psi: float, reference: str, ambient_psia: float) -> float:
    if reference == "psia":
        return value_psi
    if reference == "psig":
        return value_psi + ambient_psia
    raise ValueError("Pressure reference must be resolved to psia or psig")


def ideal_gas_properties(fluid: str, temperature_k: float) -> dict[str, float]:
    if PropsSI is None:
        raise RuntimeError("CoolProp is required for PRRV gas properties")
    gas_constant = float(PropsSI("GAS_CONSTANT", fluid) / PropsSI("MOLARMASS", fluid))
    cp = float(PropsSI("CP0MASS", "T", temperature_k, "P", 101325.0, fluid))
    return {"specific_gas_constant_j_kg_k": gas_constant, "cp_j_kg_k": cp, "cv_j_kg_k": cp - gas_constant, "gamma": cp / (cp - gas_constant)}


def ideal_gas_mass_flux(upstream_pressure_pa: float, downstream_pressure_pa: float, temperature_k: float, gas_constant: float, gamma: float) -> tuple[float, bool]:
    ratio = downstream_pressure_pa / upstream_pressure_pa
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if ratio <= critical:
        factor = math.sqrt(gamma / (gas_constant * temperature_k)) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        return upstream_pressure_pa * factor, True
    term = ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma)
    return upstream_pressure_pa * math.sqrt(2.0 * gamma * term / (gas_constant * temperature_k * (gamma - 1.0))), False


def real_gas_isentropic_mass_flux(upstream_pressure_pa: float, downstream_pressure_pa: float, upstream_temperature_k: float, fluid: str, pressure_points: int) -> dict[str, float | bool]:
    """CoolProp homogeneous-equilibrium isentropic nozzle screen."""
    if PropsSI is None:
        raise RuntimeError("CoolProp is required for high-pressure GN2 sizing")
    enthalpy_0 = float(PropsSI("H", "P", upstream_pressure_pa, "T", upstream_temperature_k, fluid))
    entropy_0 = float(PropsSI("S", "P", upstream_pressure_pa, "T", upstream_temperature_k, fluid))
    pressures = np.geomspace(downstream_pressure_pa, upstream_pressure_pa * (1.0 - 1e-10), pressure_points)
    enthalpy = np.asarray(PropsSI("H", "P", pressures, "S", entropy_0, fluid), dtype=float)
    density = np.asarray(PropsSI("D", "P", pressures, "S", entropy_0, fluid), dtype=float)
    mass_flux = density * np.sqrt(np.maximum(0.0, 2.0 * (enthalpy_0 - enthalpy)))
    index = int(np.nanargmax(mass_flux))
    return {
        "mass_flux_kg_m2_s": float(mass_flux[index]), "throat_pressure_pa": float(pressures[index]),
        "throat_to_upstream_pressure_ratio": float(pressures[index] / upstream_pressure_pa), "is_choked": bool(index > 0),
        "upstream_density_kg_m3": float(PropsSI("D", "P", upstream_pressure_pa, "T", upstream_temperature_k, fluid)),
        "upstream_compressibility_factor": float(PropsSI("Z", "P", upstream_pressure_pa, "T", upstream_temperature_k, fluid)),
    }


def calculate(config: dict[str, Any]) -> dict[str, Any]:
    """Run Mode A for known maximum flow, otherwise Mode B trade study."""
    ambient = config["ambient_pressure_psia"]
    evaluation_psia = pressure_to_psia(config["pressure"]["set_value_psi"], config["pressure"]["set_reference"], ambient)
    temperature = config["relief"]["temperature_k"]
    fluid = config["fluid"]["name"]
    real = real_gas_isentropic_mass_flux(evaluation_psia * PSI_TO_PA, config["relief"]["downstream_pressure_psia"] * PSI_TO_PA, temperature, fluid, config["analysis"]["real_gas_pressure_points"])
    mass_flux = float(real["mass_flux_kg_m2_s"])
    known_flow = config["failure"]["maximum_fill_mass_flow_kg_s"]
    if known_flow is None:
        mode = "B_fill_flow_trade_study"
        trial_flows = np.linspace(config["analysis"]["trade_fill_flow_min_kg_s"], config["analysis"]["trade_fill_flow_max_kg_s"], config["analysis"]["trade_fill_flow_points"])
    else:
        mode = "A_known_maximum_fill_flow"
        trial_flows = np.asarray([known_flow], dtype=float)
    design_flows = config["failure"]["mass_flow_factor"] * trial_flows
    required_cda = design_flows / mass_flux * M2_TO_MM2
    candidate_cda = config["candidate"]["configured_cda_mm2"]
    candidate_capacity = mass_flux * candidate_cda * MM2_TO_M2
    passes = candidate_capacity >= design_flows
    ideal = ideal_gas_properties(fluid, temperature)
    ideal_flux, ideal_choked = ideal_gas_mass_flux(evaluation_psia * PSI_TO_PA, config["relief"]["downstream_pressure_psia"] * PSI_TO_PA, temperature, ideal["specific_gas_constant_j_kg_k"], ideal["gamma"])
    # Initial COPV pressure and the ground-system flow-vs-pressure relation are
    # absent, so no constant-flow fill-time estimate is presented as physical.
    fill_time_available = False
    selected = 0 if known_flow is not None else None
    return {
        "mode": mode, "evaluation_pressure_psia": evaluation_psia, "pressure_reference_convention": config["pressure"]["set_reference"],
        "maximum_inflow_kg_s": known_flow, "design_flow_kg_s": None if selected is None else float(design_flows[selected]),
        "required_cda_mm2": None if selected is None else float(required_cda[selected]), "required_diameter_mm": None if selected is None else math.sqrt(4.0 * float(required_cda[selected]) / math.pi),
        "real_gas": real, "ideal_properties": ideal, "ideal_mass_flux_kg_m2_s": ideal_flux, "ideal_is_choked": ideal_choked,
        "candidate_cda_mm2": candidate_cda, "candidate_fully_open_envelope_kg_s": candidate_capacity,
        "trial_fill_flow_kg_s": trial_flows, "trial_design_flow_kg_s": design_flows, "trial_required_cda_mm2": required_cda,
        "trial_candidate_margin_kg_s": candidate_capacity - design_flows, "trial_candidate_pass": passes,
        "fill_time_available": fill_time_available, "fill_time_status": "UNRESOLVED_MISSING_INITIAL_STATE_AND_TRANSIENT_FILL_FLOW_MODEL",
        "estimated_fill_time_minutes": [None] * len(trial_flows),
        "candidate_analytical_status": "UNRESOLVED_MAXIMUM_FLOW" if selected is None else ("PASS" if passes[selected] else "FAIL"),
        "manufacturer_capacity_validation_status": "UNVERIFIED_NO_MANUFACTURER_RATING_OR_LIFT_CURVE",
    }


def _write_csv(path: Path, headers: list[str], rows: Iterable[Iterable[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(headers); writer.writerows(rows)


def write_outputs(config: dict[str, Any], result: dict[str, Any]) -> list[Path]:
    output = config["paths"]["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    # Remove artifacts from the superseded fixed-0.316-kg/s pressure sweep so
    # stale results cannot be mistaken for the unresolved-flow trade study.
    for legacy_name in (
        "prrv_pressure_sweep.csv",
        "prrv_property_model_audit.csv",
        "prrv_required_cda_vs_pressure.png",
    ):
        (output / legacy_name).unlink(missing_ok=True)
    written: list[Path] = []
    rows = [
        ("Requirements", "failure_scenario", config["failure"]["name"], "-", "Only LFS fill valves failed open"),
        ("Requirements", "protected_system", "two GN2 COPVs and common manifold", "-", "Ground-side filling"),
        ("Requirements", "copv_meop", config["pressure"]["meop_value_psi"], "psi", config["pressure"]["meop_reference"]),
        ("Requirements", "copv_mawp", config["pressure"]["mawp_value_psi"], "psi", config["pressure"]["mawp_reference"]),
        ("Requirements", "prrv_cracking_pressure_target", config["pressure"]["set_value_psi"], "psi", "Derived from configured PRRV crack pressure and ambient; conflicts with documented MAWP"),
        ("Requirements", "evaluated_relieving_pressure", result["evaluation_pressure_psia"], "psia", "Setpoint converted explicitly"),
        ("Inputs", "maximum_LFS_fill_flow", result["maximum_inflow_kg_s"], "kg/s", config["failure"]["maximum_fill_flow_status"]),
        ("Inputs", "analysis_mode", result["mode"], "-", "Known-flow mode is automatic when maximum flow is supplied"),
        ("Inputs", "gn2_relieving_temperature", config["relief"]["temperature_k"], "K", "Minimum configured COPV temperature"),
        ("Inputs", "atmospheric_discharge_pressure", config["relief"]["downstream_pressure_psia"], "psia", "hopper_config.json"),
        ("Inputs", "total_copv_volume", sum(config["copvs"]["volume_l"]), "L", "Two repository COPVs"),
        ("Inputs", "candidate_valve_cda", result["candidate_cda_mm2"], "mm^2", "Configured guess; not certified"),
        ("Inputs", "PRQD_or_LFS_orifice_cda", None, "mm^2", "Unresolved; not invented"),
        ("Assumptions", "mass_flow_factor", config["failure"]["mass_flow_factor"], "-", "Applied exactly once to each trial flow"),
        ("Assumptions", "gas_property_model", "CoolProp real-gas Nitrogen", "-", "293 K"),
        ("Assumptions", "steady_state_limitations", "Fully-open theoretical capacity only", "-", "No immediate full capacity assumed at crack"),
        ("Methodology", "required_cda", f"{config['failure']['mass_flow_factor']:.6g} * trial LFS flow / real-gas mass flux", "-", "Choked state determined along isentrope"),
        ("Results", "maximum_LFS_fill_flow_status", config["failure"]["maximum_fill_flow_status"], "-", "0.316 kg/s is not treated as authoritative"),
        ("Results", "required_design_relief_flow", result["design_flow_kg_s"], "kg/s", "Blank until maximum flow is supplied"),
        ("Results", "required_effective_cda", result["required_cda_mm2"], "mm^2", "Blank in unresolved-flow Mode B"),
        ("Results", "candidate_theoretical_mass_flow_capacity", result["candidate_fully_open_envelope_kg_s"], "kg/s", f"At {result['evaluation_pressure_psia']:.1f} psia and {config['relief']['temperature_k']:.1f} K"),
        ("Results", "candidate_screening", result["candidate_analytical_status"], "-", "Per-flow results in trade CSV"),
        ("Results", "fill_time_status", result["fill_time_status"], "-", "No volume-only or unjustified constant-flow estimate"),
        ("Results", "final_validation_status", "UNRESOLVED", "-", "LFS flow, MAWP conflict, and certified valve data unresolved"),
        ("Results", "missing_inputs_preventing_final_sizing", "; ".join(config["missing_inputs"]), "-", "Engineering action list"),
    ]
    summary = output / "prrv_sizing_summary.csv"
    _write_csv(summary, ["section", "metric", "value", "unit", "note"], rows); written.append(summary)
    trade = output / "prrv_fill_flow_trade_study.csv"
    _write_csv(trade, ["assumed_LFS_fill_flow_kg_s", "design_relief_flow_kg_s", "required_effective_cda_mm2", "candidate_theoretical_capacity_kg_s", "candidate_margin_kg_s", "candidate_pass", "estimated_fill_time_minutes", "meets_30_min_fill_target"], zip(result["trial_fill_flow_kg_s"], result["trial_design_flow_kg_s"], result["trial_required_cda_mm2"], np.full_like(result["trial_fill_flow_kg_s"], result["candidate_fully_open_envelope_kg_s"]), result["trial_candidate_margin_kg_s"], result["trial_candidate_pass"].astype(int), result["estimated_fill_time_minutes"], np.zeros_like(result["trial_fill_flow_kg_s"], dtype=int))); written.append(trade)
    unresolved = output / "prrv_unresolved_inputs.csv"
    _write_csv(unresolved, ["unresolved_engineering_input"], ((item,) for item in config["missing_inputs"])); written.append(unresolved)
    conflicts = output / "prrv_configuration_conflicts.csv"
    _write_csv(conflicts, ["configuration_conflict"], ((item,) for item in config["configuration_conflicts"])); written.append(conflicts)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(result["trial_fill_flow_kg_s"], result["trial_required_cda_mm2"], linewidth=2, label="Required effective CdA")
    axis.axhline(result["candidate_cda_mm2"], color="tab:orange", linestyle="--", label="Configured candidate CdA")
    axis.set(xlabel="Assumed maximum LFS fill flow [kg/s]", ylabel="Required effective CdA [mm²]", title="PRRV unresolved-fill-flow trade study")
    axis.grid(True); axis.legend(); path = output / "prrv_required_cda_vs_fill_flow.png"; fig.savefig(path, dpi=200); plt.close(fig); written.append(path)

    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    axis.plot(result["trial_fill_flow_kg_s"], result["trial_design_flow_kg_s"], linewidth=2, label=f"Required design flow ({config['failure']['mass_flow_factor']:.3g}×)")
    axis.axhline(result["candidate_fully_open_envelope_kg_s"], color="tab:orange", linestyle="--", label="Candidate fully-open theoretical capacity")
    axis.set(xlabel="Assumed maximum LFS fill flow [kg/s]", ylabel="GN2 mass flow [kg/s]", title="PRRV candidate capacity trade — not certified")
    axis.grid(True); axis.legend(); path = output / "prrv_candidate_capacity_comparison.png"; fig.savefig(path, dpi=200); plt.close(fig); written.append(path)
    return written


def print_report(config: dict[str, Any], result: dict[str, Any], written: list[Path]) -> None:
    print("\nPRELIMINARY PRRV LFS FILL-VALVES-FAILED-OPEN SCREEN")
    print(f"Mode: {result['mode']}; maximum LFS flow: {result['maximum_inflow_kg_s']}")
    print(f"{result['evaluation_pressure_psia']:.1f} psia candidate theoretical capacity: {result['candidate_fully_open_envelope_kg_s']:.6f} kg/s")
    print(f"Trade range: {result['trial_fill_flow_kg_s'][0]:.3f} to {result['trial_fill_flow_kg_s'][-1]:.3f} kg/s")
    print(f"Fill time: {result['fill_time_status']}")
    for item in config["configuration_conflicts"]: print(f"WARNING: {item}")
    for path in written: print(f"  - {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_HOPPER_CONFIG)
    config = load_prrv_inputs(parser.parse_args().config).as_solver_config()
    warnings.warn("PRRV results are theoretical screening, not certified relief capacity", stacklevel=2)
    result = calculate(config); written = write_outputs(config, result); print_report(config, result, written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

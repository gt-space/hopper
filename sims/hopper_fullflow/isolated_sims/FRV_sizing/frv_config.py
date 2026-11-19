"""Load and validate inputs for the preliminary FRV analytical sizer.

Every vehicle, hardware, and sizing input is read directly from
``hopper_config.json`` on every call. This module contains no sizing equations
and never writes the input file.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from CoolProp.CoolProp import PropsSI


HERE = Path(__file__).resolve().parent
DEFAULT_HOPPER_CONFIG = HERE.parents[1] / "hopper_config.json"


class FRVConfigError(ValueError):
    """An FRV input is absent, ambiguous, or physically invalid."""


# Public attribute name -> (validated solver-config path, hopper_config.json key).
# Values are never stored here; this is only the traceable wiring diagram.
FRV_INPUT_SOURCES: dict[str, tuple[tuple[str, ...], str]] = {
    "hopper_config_path": (("paths", "hopper_config"), "hopper_config.json"),
    "output_directory": (("paths", "output_directory"), "files.results_directory + relief_sizing.frv.output_subdirectory"),
    "pressurant_name": (("pressurant", "name"), "fluids.pressurant"),
    "pressurant_model": (("pressurant", "model"), "relief_sizing.frv.pressurant_model"),
    "specific_gas_constant_j_kg_k": (("pressurant", "specific_gas_constant_j_kg_k"), "CoolProp(fluids.pressurant, bottles[*])"),
    "specific_heat_ratio": (("pressurant", "specific_heat_ratio"), "CoolProp(fluids.pressurant, bottles[*])"),
    "ambient_pressure_psia": (("ambient_pressure_psia",), "fluids.ambient_pressure_psia"),
    "source_names": (("source", "names"), "bottles[*].name"),
    "source_volumes_l": (("source", "volume_l"), "bottles[*].volume_l"),
    "source_pressures_psia": (("source", "pressure_psia"), "bottles[*].pressure_psia"),
    "source_temperatures_k": (("source", "temperature_k"), "bottles[*].temperature_k"),
    "source_outlet_cdas_mm2": (("source", "outlet_cda_mm2"), "bottles[*].outlet_cda_mm2"),
    "source_model": (("source", "model"), "pressurization.copv_model"),
    "source_header_volume_l": (("source", "header_volume_l"), "pressurization.header_volume_l"),
    "tank_name": (("tank", "name"), "tanks.fuel.name"),
    "tank_propellant": (("tank", "propellant"), "tanks.fuel.propellant"),
    "tank_total_volume_l": (("tank", "total_volume_l"), "tanks.fuel.volume_l"),
    "tank_initial_liquid_mass_kg": (("tank", "initial_liquid_mass_kg"), "tanks.fuel.liquid_mass_kg"),
    "tank_initial_pressure_psia": (("tank", "initial_pressure_psia"), "tanks.fuel.pressure_psia"),
    "tank_ullage_temperature_k": (("tank", "ullage_temperature_k"), "tanks.fuel.ullage_temperature_k"),
    "tank_control_set_pressure_psia": (("tank", "control_set_pressure_psia"), "tanks.fuel.set_pressure_psia"),
    "tank_meop_psia": (("tank", "meop_psia"), "tanks.fuel.set_pressure_psia"),
    "tank_mawp_psia": (("tank", "mawp_psia"), "relief_sizing.frv.tank_mawp_psia"),
    "maximum_relieving_pressure_psia": (("tank", "maximum_allowable_relieving_pressure_psia"), "side_paths[FRV Relief].elements[FRV].full_open_pressure_psig + fluids.ambient_pressure_psia"),
    "failure_scenario": (("failure", "name"), "relief_sizing.frv.failure_scenario"),
    "simulation_maximum_gn2_inflow_kg_s": (("failure", "simulation_maximum_gn2_inflow_kg_s"), "relief_sizing.frv.simulation_maximum_gn2_inflow_kg_s"),
    "simulation_inflow_status": (("failure", "simulation_inflow_status"), "relief_sizing.frv.simulation_inflow_status"),
    "pressurant_line": (("failure", "pressurant_line"), "flow_paths.fuel_pressurant[Fuel Pressurant Line]"),
    "fbb_cda_mm2": (("failure", "fbb_cda_mm2"), "flow_paths.fuel_pressurant[FBB].cda_mm2"),
    "fbb_open_time_s": (("failure", "fbb_open_time_s"), "flow_paths.fuel_pressurant[FBB].open_time_s"),
    "fbb_close_time_s": (("failure", "fbb_close_time_s"), "flow_paths.fuel_pressurant[FBB].close_time_s"),
    "orifice_name": (("failure", "orifice_name"), "flow_paths.fuel_pressurant[OR-2 or FBBORF].name"),
    "orifice_cda_mm2": (("failure", "orifice_cda_mm2"), "flow_paths.fuel_pressurant[OR-2 or FBBORF].cda_mm2"),
    "check_valve_cda_mm2": (("failure", "check_valve_cda_mm2"), "flow_paths.fuel_pressurant[PRFCK].cda_mm2"),
    "sizing_margin_fraction": (("failure", "sizing_margin_fraction"), "relief_sizing.shared.mass_flow_safety_factor - 1"),
    "frv_name": (("frv", "name"), "side_paths[FRV Relief].elements[FRV].name"),
    "frv_source": (("frv", "source"), "side_paths[FRV Relief].source"),
    "frv_set_pressure_psia": (("frv", "proposed_set_pressure_psia"), "side_paths[FRV Relief].elements[FRV].crack_pressure_psig + fluids.ambient_pressure_psia"),
    "frv_crack_pressure_psig": (("frv", "configured_hardware_crack_pressure_psig"), "side_paths[FRV Relief].elements[FRV].crack_pressure_psig"),
    "frv_full_open_cda_mm2": (("frv", "configured_full_open_cda_mm2"), "side_paths[FRV Relief].elements[FRV].cda_mm2"),
    "frv_lift_model": (("frv", "lift_model"), "relief_sizing.frv.lift_model"),
    "regulator_name": (("context", "regulator_name"), "regulators[PR-REG or PUREG].name"),
    "regulator_set_pressure_psig": (("context", "regulator_set_pressure_psig"), "regulators[PR-REG or PUREG].set_pressure_psig"),
    "regulator_cda_mm2": (("context", "regulator_cda_mm2"), "regulators[PR-REG or PUREG].cda_mm2"),
    "analysis_mode": (("analysis", "mode"), "relief_sizing.frv.analysis_mode"),
    "pressure_sweep_points": (("analysis", "pressure_sweep_points"), "relief_sizing.shared.pressure_sweep_points"),
}


@dataclass(frozen=True)
class FRVInputs:
    """Named, read-only view of one freshly loaded JSON configuration."""

    _solver_config: dict[str, Any]

    def __getattr__(self, name: str) -> Any:
        try:
            path, _ = FRV_INPUT_SOURCES[name]
        except KeyError as exc:
            raise AttributeError(name) from exc
        value: Any = self._solver_config
        for key in path:
            value = value[key]
        return value

    def as_solver_config(self) -> dict[str, Any]:
        """Return the validated mapping used by the unchanged calculations."""
        return self._solver_config


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise FRVConfigError(f"Duplicate JSON field: {key}")
        value[key] = item
    return value


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise FRVConfigError(f"Could not find {label}: {path}")
    try:
        value = json.loads(
            path.read_text(encoding="utf-8-sig"),
            object_pairs_hook=_unique_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                FRVConfigError(f"{label} contains invalid number {token}")
            ),
        )
    except json.JSONDecodeError as exc:
        raise FRVConfigError(f"Could not parse {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise FRVConfigError(f"{label} must contain a JSON object at its root")
    return value


def _require_fields(value: Any, fields: tuple[str, ...], label: str) -> None:
    if not isinstance(value, dict):
        raise FRVConfigError(f"{label} must be a JSON object")
    missing = [field for field in fields if field not in value]
    empty = [field for field in fields if field in value and value[field] is None]
    if missing or empty:
        raise FRVConfigError(
            f"{label}: missing fields {missing}; null required fields {empty}"
        )


def _named(items: Any, wanted_name: str, label: str) -> dict[str, Any]:
    if not isinstance(items, list) or not all(isinstance(item, dict) for item in items):
        raise FRVConfigError(f"{label} must be an array of JSON objects")
    matches = [item for item in items if item.get("name") == wanted_name]
    if len(matches) != 1:
        raise FRVConfigError(
            f'{label} must contain exactly one object named "{wanted_name}"; '
            f"found {len(matches)}"
        )
    return matches[0]


def _named_any(items: Any, wanted_names: tuple[str, ...], label: str) -> dict[str, Any]:
    """Resolve one component across approved P&ID renames, rejecting ambiguity."""
    if not isinstance(items, list) or not all(isinstance(item, dict) for item in items):
        raise FRVConfigError(f"{label} must be an array of JSON objects")
    matches = [item for item in items if item.get("name") in wanted_names]
    if len(matches) != 1:
        raise FRVConfigError(
            f"{label} must contain exactly one of {list(wanted_names)}; found {len(matches)}"
        )
    return matches[0]


def _nitrogen_properties(fluid: str, temperature_k: float, pressure_psia: float) -> tuple[float, float]:
    """Derive ideal-gas R and gamma from the configured fluid and source state."""
    pressure_pa = pressure_psia * 6894.757293168
    gas_constant = float(PropsSI("GAS_CONSTANT", fluid) / PropsSI("MOLAR_MASS", fluid))
    cp = float(PropsSI("CP0MASS", "T", temperature_k, "P", pressure_pa, fluid))
    return gas_constant, cp / (cp - gas_constant)


def _nullable_number(value: Any, label: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise FRVConfigError(f"{label} must be a finite number or null")
    if not math.isfinite(value):
        raise FRVConfigError(f"{label} must be finite or null")
    return float(value)


def _input_mapping() -> list[dict[str, str]]:
    """Trace every Python input to the exact source field and source units."""
    rows = [
        ("pressurant.name", "fluids.pressurant", "-", "Pressurant identity"),
        ("ambient_pressure_psia", "fluids.ambient_pressure_psia", "psia", "FRV outlet pressure and gauge comparison"),
        ("source.model", "pressurization.copv_model", "-", "COPV topology"),
        ("source.header_volume_l", "pressurization.header_volume_l", "L", "Shared-header context; unused in lumped steady case"),
        ("source.volume_l", "bottles[*].volume_l", "L", "COPV inventory context"),
        ("source.pressure_psia", "bottles[*].pressure_psia", "psia", "FBB upstream stagnation pressure"),
        ("source.temperature_k", "bottles[*].temperature_k", "K", "FBB upstream stagnation temperature"),
        ("source.outlet_cda_mm2", "bottles[*].outlet_cda_mm2", "mm^2", "Parallel COPV outlet restrictions"),
        ("tank.initial_pressure_psia", "tanks.fuel.pressure_psia", "psia", "Initial tank-pressure check"),
        ("tank.control_set_pressure_psia", "tanks.fuel.set_pressure_psia", "psia", "Bang-bang control setpoint"),
        ("tank.meop_psia", "tanks.fuel.set_pressure_psia", "psia", "Authoritative 500 psia MEOP source per engineering direction"),
        ("tank.total_volume_l", "tanks.fuel.volume_l", "L", "Future transient context"),
        ("tank.initial_liquid_mass_kg", "tanks.fuel.liquid_mass_kg", "kg", "Future ullage-volume context"),
        ("tank.ullage_temperature_k", "tanks.fuel.ullage_temperature_k", "K", "Preliminary FRV upstream temperature"),
        ("failure.pressurant_line.length_in", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].length_in", "in", "Pressurant-line length"),
        ("failure.pressurant_line.inner_diameter_in", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].inner_diameter_in", "in", "Optional direct ID; null means derive from OD and wall"),
        ("failure.pressurant_line.outer_diameter_in", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].outer_diameter_in", "in", "Pressurant-line OD"),
        ("failure.pressurant_line.wall_thickness_in", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].wall_thickness_in", "in", "Pressurant-line wall thickness"),
        ("failure.pressurant_line.roughness_in", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].roughness_in", "in", "Fully rough line-CdA approximation"),
        ("failure.pressurant_line.inertia", "flow_paths.fuel_pressurant[name=Fuel Pressurant Line].inertia", "boolean", "Selects steady equivalent-CdA line model"),
        ("failure.fbb_cda_mm2", "flow_paths.fuel_pressurant[name=FBB].cda_mm2", "mm^2", "Failed-open valve restriction"),
        ("failure.fbb_open_time_s", "flow_paths.fuel_pressurant[name=FBB].open_time_s", "s", "Future opening transient"),
        ("failure.fbb_close_time_s", "flow_paths.fuel_pressurant[name=FBB].close_time_s", "s", "Future valve-dynamics context"),
        ("failure.orifice_cda_mm2", "flow_paths.fuel_pressurant[name=OR-2 or FBBORF].cda_mm2", "mm^2", "Downstream orifice restriction"),
        ("failure.check_valve_cda_mm2", "flow_paths.fuel_pressurant[name=PRFCK].cda_mm2", "mm^2", "Downstream check-valve restriction"),
        ("failure.simulation_maximum_gn2_inflow_kg_s", "relief_sizing.frv.simulation_maximum_gn2_inflow_kg_s", "kg/s", "Authoritative latest free-flight maximum"),
        ("frv.configured_hardware_crack_pressure_psig", "side_paths[name=FRV Relief].elements[name=FRV].crack_pressure_psig", "psig", "Existing hardware comparison only"),
        ("frv.configured_full_open_cda_mm2", "side_paths[name=FRV Relief].elements[name=FRV].cda_mm2", "mm^2", "Analytical fully-open envelope only"),
        ("context.regulator_set_pressure_psig", "regulators[name=PR-REG or PUREG].set_pressure_psig", "psig", "Purge-regulator context; excluded from FBB path"),
        ("context.regulator_cda_mm2", "regulators[name=PR-REG or PUREG].cda_mm2", "mm^2", "Purge-regulator context; excluded from FBB path"),
        ("frv.proposed_set_pressure_psia", "FRV crack_pressure_psig + fluids.ambient_pressure_psia", "psia", "Configured relief set pressure"),
        ("tank.mawp_psia", "relief_sizing.frv.tank_mawp_psia", "psia", "Tank MAWP requirement"),
        ("tank.maximum_allowable_relieving_pressure_psia", "FRV full_open_pressure_psig + ambient (default 1.1 * crack gauge)", "psia", "Configured full-open pressure"),
        ("failure.name", "relief_sizing.frv.failure_scenario", "-", "Primary analytical failure case"),
        ("failure.sizing_margin_fraction", "relief_sizing.shared.mass_flow_safety_factor - 1", "-", "Positive sizing margin M"),
        ("pressurant.properties", "CoolProp from fluids.pressurant and configured COPV state", "SI", "Ideal-gas equation properties"),
        ("analysis.pressure_sweep_points", "relief_sizing.shared.pressure_sweep_points", "-", "Pressure-grid resolution"),
        ("analysis.mode", "relief_sizing.frv.analysis_mode", "-", "Steady-state-only analysis selection"),
        ("paths.output_directory", "files.results_directory + relief_sizing.frv.output_subdirectory", "-", "FRV CSV and plot destination"),
    ]
    return [
        {"python_variable": variable, "json_source_path": source, "source_units": units, "purpose": purpose}
        for variable, source, units, purpose in rows
    ]


def load_frv_config(
    config_path: str | Path = DEFAULT_HOPPER_CONFIG,
) -> dict[str, Any]:
    """Return validated, traceable FRV inputs without modifying source JSON."""
    config_path = Path(config_path).expanduser().resolve()
    raw = _read_json(config_path, "Hopper configuration")

    _require_fields(raw, ("files", "fluids", "pressurization", "bottles", "tanks", "flow_paths", "regulators", "side_paths", "relief_sizing"), "hopper_config.json")
    _require_fields(raw["files"], ("results_directory",), "files")
    _require_fields(raw["fluids"], ("pressurant", "ambient_pressure_psia"), "fluids")
    _require_fields(raw["pressurization"], ("copv_model", "header_volume_l"), "pressurization")
    _require_fields(raw["tanks"], ("fuel",), "tanks")
    _require_fields(raw["flow_paths"], ("fuel_pressurant",), "flow_paths")

    _require_fields(raw["relief_sizing"], ("shared", "frv"), "relief_sizing")
    shared = raw["relief_sizing"]["shared"]
    assumptions = raw["relief_sizing"]["frv"]
    _require_fields(shared, ("mass_flow_safety_factor", "pressure_sweep_points"), "relief_sizing.shared")
    _require_fields(assumptions, ("failure_scenario", "simulation_maximum_gn2_inflow_kg_s", "simulation_inflow_status", "tank_mawp_psia", "pressurant_model", "lift_model", "analysis_mode", "output_subdirectory"), "relief_sizing.frv")

    fuel_tank = raw["tanks"]["fuel"]
    fuel_path = raw["flow_paths"]["fuel_pressurant"]
    pressurant_line = _named(fuel_path, "Fuel Pressurant Line", "flow_paths.fuel_pressurant")
    fbb = _named(fuel_path, "FBB", "flow_paths.fuel_pressurant")
    orifice = _named_any(fuel_path, ("OR-2", "FBBORF"), "flow_paths.fuel_pressurant")
    check_valve = _named(fuel_path, "PRFCK", "flow_paths.fuel_pressurant")
    frv_branch = _named(raw["side_paths"], "FRV Relief", "side_paths")
    _require_fields(frv_branch, ("source", "elements"), "side_paths[name=FRV Relief]")
    frv = _named(frv_branch["elements"], "FRV", "side_paths[name=FRV Relief].elements")
    regulator = _named_any(raw["regulators"], ("PR-REG", "PUREG"), "regulators")

    _require_fields(fuel_tank, ("name", "propellant", "volume_l", "liquid_mass_kg", "pressure_psia", "ullage_temperature_k", "set_pressure_psia"), "tanks.fuel")
    _require_fields(pressurant_line, ("length_in", "roughness_in", "inertia"), "Fuel Pressurant Line")
    if pressurant_line.get("inner_diameter_in") is None:
        _require_fields(pressurant_line, ("outer_diameter_in", "wall_thickness_in"), "Fuel Pressurant Line")
    _require_fields(fbb, ("cda_mm2", "open_time_s", "close_time_s"), "FBB")
    _require_fields(orifice, ("cda_mm2",), orifice["name"])
    _require_fields(check_valve, ("cda_mm2",), "PRFCK")
    _require_fields(frv, ("cda_mm2", "crack_pressure_psig"), "FRV")
    _require_fields(regulator, ("set_pressure_psig", "cda_mm2"), regulator["name"])

    bottles = raw["bottles"]
    if not isinstance(bottles, list) or not bottles:
        raise FRVConfigError("hopper_config.json must define at least one COPV")
    for index, bottle in enumerate(bottles):
        _require_fields(bottle, ("name", "volume_l", "pressure_psia", "temperature_k", "outlet_cda_mm2"), f"bottles[{index}]")

    ambient = raw["fluids"]["ambient_pressure_psia"]
    crack_psig = frv["crack_pressure_psig"]
    full_open_psig = frv.get("full_open_pressure_psig") or 1.1 * crack_psig
    set_pressure_psia = crack_psig + ambient
    maximum_relieving_pressure_psia = full_open_psig + ambient
    gas_constant, gamma = _nitrogen_properties(raw["fluids"]["pressurant"], min(b["temperature_k"] for b in bottles), max(b["pressure_psia"] for b in bottles))
    output_directory = (
        config_path.parent / raw["files"]["results_directory"] / assumptions["output_subdirectory"]
    ).resolve()
    config: dict[str, Any] = {
        "paths": {
            "hopper_config": config_path,
            "output_directory": output_directory,
        },
        "source_paths": {
            "tank_meop_psia": "tanks.fuel.set_pressure_psia",
            "tank_initial_pressure_psia": "tanks.fuel.pressure_psia",
        },
        "pressurant": {
            "name": raw["fluids"]["pressurant"],
            "model": assumptions["pressurant_model"],
            "specific_gas_constant_j_kg_k": gas_constant,
            "specific_heat_ratio": gamma,
        },
        "ambient_pressure_psia": ambient,
        "source": {
            "names": [bottle["name"] for bottle in bottles],
            "volume_l": [bottle["volume_l"] for bottle in bottles],
            "pressure_psia": [bottle["pressure_psia"] for bottle in bottles],
            "temperature_k": [bottle["temperature_k"] for bottle in bottles],
            "outlet_cda_mm2": [bottle["outlet_cda_mm2"] for bottle in bottles],
            "model": raw["pressurization"]["copv_model"],
            "header_volume_l": raw["pressurization"]["header_volume_l"],
        },
        "tank": {
            "name": fuel_tank["name"],
            "propellant": fuel_tank["propellant"],
            "total_volume_l": fuel_tank["volume_l"],
            "initial_liquid_mass_kg": fuel_tank["liquid_mass_kg"],
            "initial_pressure_psia": fuel_tank["pressure_psia"],
            "ullage_temperature_k": fuel_tank["ullage_temperature_k"],
            "control_set_pressure_psia": fuel_tank["set_pressure_psia"],
            "meop_psia": fuel_tank["set_pressure_psia"],
            "mawp_psia": _nullable_number(assumptions["tank_mawp_psia"], "tank_mawp_psia"),
            "mawp_status": "relief_sizing.frv.tank_mawp_psia",
            "permitted_pressure_accumulation_psi": maximum_relieving_pressure_psia - set_pressure_psia,
            "maximum_allowable_relieving_pressure_psia": maximum_relieving_pressure_psia,
            "relieving_pressure_status": "derived_from_configured_relief_full_open_pressure",
        },
        "failure": {
            "name": assumptions["failure_scenario"],
            "simulation_maximum_gn2_inflow_kg_s": assumptions["simulation_maximum_gn2_inflow_kg_s"],
            "simulation_inflow_status": assumptions["simulation_inflow_status"],
            "pressurant_line": dict(pressurant_line),
            "fbb_cda_mm2": fbb["cda_mm2"],
            "fbb_open_time_s": fbb["open_time_s"],
            "fbb_close_time_s": fbb["close_time_s"],
            "orifice_cda_mm2": orifice["cda_mm2"],
            "orifice_name": orifice["name"],
            "check_valve_cda_mm2": check_valve["cda_mm2"],
            "sizing_margin_fraction": shared["mass_flow_safety_factor"] - 1.0,
            "simultaneous_fbb_obb_requires_separate_transient_evaluation": True,
        },
        "frv": {
            "name": frv["name"],
            "source": frv_branch["source"],
            "proposed_set_pressure_psia": set_pressure_psia,
            "set_pressure_status": "derived_from_configured_hardware_crack_pressure",
            "configured_hardware_crack_pressure_psig": frv["crack_pressure_psig"],
            "configured_full_open_cda_mm2": frv["cda_mm2"],
            "lift_model": assumptions["lift_model"],
            "full_lift_overpressure_fraction": (full_open_psig - crack_psig) / crack_psig,
            "certified_capacity_available": False,
        },
        "context": {
            "regulator_name": regulator["name"],
            "regulator_set_pressure_psig": regulator["set_pressure_psig"],
            "regulator_cda_mm2": regulator["cda_mm2"],
            "regulator_used_in_fbb_failure_path": False,
        },
        "analysis": {
            "mode": assumptions["analysis_mode"],
            "pressure_sweep_points": shared["pressure_sweep_points"],
            "pressure_sweep_upper_basis": "configured_relief_full_open_pressure",
            "transient_enabled": False,
            "preliminary_only": True,
        },
        "input_mapping": _input_mapping(),
        "missing_inputs": [
            "Lead/code approval of relief_sizing.frv.tank_mawp_psia",
            "Lead/code approval of the configured FRV crack and full-open pressures",
            "Manufacturer FRV full-lift pressure or certified capacity curve",
            "Manufacturer-certified GN2 relieving capacity and applicable coefficient",
            "FRV inlet/outlet plumbing geometry and maximum backpressure",
            "Validated real-gas relief method/properties at the COPV condition",
            "Validated multirestriction compressible-flow method or intermediate pressures",
            "Validated flow data for COPV outlets, FBB, FBBORF/OR-2, and PRFCK",
            "Fuel-ullage relieving temperature",
            "Fuel ullage volume, thermal model, valve dynamics, timestep, and end time for a transient",
        ],
    }
    validate_frv_config(config)
    return config


def load_frv_inputs(config_path: str | Path = DEFAULT_HOPPER_CONFIG) -> FRVInputs:
    """Reload JSON and expose every FRV input through a named attribute."""
    return FRVInputs(load_frv_config(config_path))


def validate_frv_config(config: dict[str, Any]) -> None:
    """Validate types, units, topology, and supported preliminary models."""
    source = config["source"]
    tank = config["tank"]
    failure = config["failure"]
    frv = config["frv"]
    pressurant = config["pressurant"]
    analysis = config["analysis"]

    positive = [
        config["ambient_pressure_psia"], *source["volume_l"], *source["pressure_psia"],
        *source["temperature_k"], *source["outlet_cda_mm2"], tank["initial_pressure_psia"],
        tank["control_set_pressure_psia"], tank["meop_psia"], tank["ullage_temperature_k"],
        failure["simulation_maximum_gn2_inflow_kg_s"], failure["fbb_cda_mm2"], failure["orifice_cda_mm2"], failure["check_valve_cda_mm2"],
        failure["sizing_margin_fraction"], frv["proposed_set_pressure_psia"],
        frv["configured_hardware_crack_pressure_psig"], frv["configured_full_open_cda_mm2"],
        frv["full_lift_overpressure_fraction"], pressurant["specific_gas_constant_j_kg_k"],
        pressurant["specific_heat_ratio"], analysis["pressure_sweep_points"],
    ]
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 for value in positive):
        raise FRVConfigError("Required FRV numerical inputs must be finite and positive")
    if pressurant["specific_heat_ratio"] <= 1:
        raise FRVConfigError("Specific heat ratio must be greater than one")
    if failure["name"] != "FBB_fail_open_alone":
        raise FRVConfigError(f"Unsupported failure scenario: {failure['name']}")
    if failure["sizing_margin_fraction"] <= 0:
        raise FRVConfigError("The configured FRV sizing margin must be positive")
    if analysis["mode"] != "steady_state" or analysis["transient_enabled"]:
        raise FRVConfigError("This solver implements steady-state analysis only")
    if analysis["pressure_sweep_upper_basis"] != "configured_relief_full_open_pressure":
        raise FRVConfigError("Unsupported pressure-sweep upper basis")
    if not isinstance(analysis["pressure_sweep_points"], int) or isinstance(analysis["pressure_sweep_points"], bool) or analysis["pressure_sweep_points"] < 2:
        raise FRVConfigError("pressure_sweep_points must be an integer of at least two")
    if pressurant["model"] != "ideal_gas_constant_gamma":
        raise FRVConfigError(f"Unsupported pressurant model: {pressurant['model']}")
    if frv["lift_model"] != "configured_relief_valve_linear":
        raise FRVConfigError(f"Unsupported FRV lift model: {frv['lift_model']}")

    line = failure["pressurant_line"]
    if line["inertia"] is not False:
        raise FRVConfigError("Steady sizing requires the configured non-inertial equivalent-CdA line model")
    diameter_in = line.get("inner_diameter_in")
    if diameter_in is None:
        diameter_in = line["outer_diameter_in"] - 2 * line["wall_thickness_in"]
    for value, label in ((diameter_in, "inside diameter"), (line["length_in"], "length"), (line["roughness_in"], "roughness")):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise FRVConfigError(f"Fuel pressurant line {label} must be finite and positive")
    if line["roughness_in"] >= diameter_in:
        raise FRVConfigError("Fuel pressurant line roughness must be less than its inside diameter")

    if source["model"] == "lumped":
        if max(source["pressure_psia"]) - min(source["pressure_psia"]) > 1e-9 or max(source["temperature_k"]) - min(source["temperature_k"]) > 1e-9:
            raise FRVConfigError("Lumped COPVs must have equal initial pressure and temperature")
    if tank["maximum_allowable_relieving_pressure_psia"] is not None:
        if tank["maximum_allowable_relieving_pressure_psia"] <= frv["proposed_set_pressure_psia"]:
            raise FRVConfigError("Maximum relieving pressure must exceed proposed set pressure")
        if tank["mawp_psia"] is not None and tank["maximum_allowable_relieving_pressure_psia"] > tank["mawp_psia"]:
            raise FRVConfigError("Maximum relieving pressure exceeds configured MAWP")
        accumulation = tank["permitted_pressure_accumulation_psi"]
        if accumulation is None or accumulation <= 0:
            raise FRVConfigError("Known maximum relieving pressure requires positive permitted accumulation")
        expected_relief_pressure = frv["proposed_set_pressure_psia"] + accumulation
        if not math.isclose(
            tank["maximum_allowable_relieving_pressure_psia"], expected_relief_pressure,
            rel_tol=0.0, abs_tol=1e-9,
        ):
            raise FRVConfigError(
                "Maximum relieving pressure must equal set pressure plus permitted accumulation"
            )


if __name__ == "__main__":
    loaded = load_frv_config()
    print(f"Validated Hopper inputs: {loaded['paths']['hopper_config']}")
    print("Validated FRV sizing inputs: hopper_config.json:relief_sizing.frv")

"""Input boundary for preliminary oxidizer B/B fail-open ORV sizing.

Hardware and operating-state values are read from ``hopper_config.json``. This
module contains only ORV-specific requirements or assumptions absent there.
The calculation is steady-state and gas-only; transient inputs are not loaded.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
FULLFLOW_ROOT = HERE.parents[1]
DEFAULT_HOPPER_CONFIG = FULLFLOW_ROOT / "hopper_config.json"

class ORVConfigError(ValueError):
    """An ORV input is absent, ambiguous, or physically inconsistent."""


ORV_INPUT_SOURCES: dict[str, tuple[tuple[str, ...], str]] = {
    "hopper_config_path": (("paths", "hopper_config"), "hopper_config.json"),
    "output_directory": (("paths", "output_directory"), "files.results_directory + relief_sizing.orv.output_subdirectory"),
    "pressurant_name": (("pressurant", "name"), "fluids.pressurant"),
    "ambient_pressure_psia": (("ambient_pressure_psia",), "fluids.ambient_pressure_psia"),
    "source_names": (("source", "names"), "bottles[*].name"),
    "source_pressures_psia": (("source", "pressure_psia"), "bottles[*].pressure_psia"),
    "source_temperatures_k": (("source", "temperature_k"), "bottles[*].temperature_k"),
    "source_outlet_cdas_mm2": (("source", "outlet_cda_mm2"), "bottles[*].outlet_cda_mm2"),
    "tank_name": (("tank", "name"), "tanks.oxidizer.name"),
    "tank_propellant": (("tank", "propellant"), "tanks.oxidizer.propellant"),
    "tank_initial_pressure_psia": (("tank", "initial_pressure_psia"), "tanks.oxidizer.pressure_psia"),
    "tank_liquid_temperature_k": (("tank", "liquid_temperature_k"), "tanks.oxidizer.liquid_temperature_k"),
    "tank_meop_psia": (("tank", "meop_psia"), "tanks.oxidizer.set_pressure_psia"),
    "tank_mawp_psia": (("tank", "mawp_psia"), "relief_sizing.orv.tank_mawp_psia"),
    "provisional_relieving_temperature_k": (("thermal", "provisional_relieving_temperature_k"), "tanks.oxidizer.ullage_temperature_k"),
    "liquid_to_ullage_ua_w_k": (("thermal", "liquid_to_ullage_ua_w_k"), "relief_sizing.orv.liquid_to_ullage_ua_w_k"),
    "wall_to_ullage_ua_w_k": (("thermal", "wall_to_ullage_ua_w_k"), "relief_sizing.orv.wall_to_ullage_ua_w_k"),
    "tank_wall_temperature_k": (("thermal", "tank_wall_temperature_k"), "relief_sizing.orv.tank_wall_temperature_k"),
    "net_external_sensible_heat_rate_w": (("thermal", "net_external_sensible_heat_rate_w"), "relief_sizing.orv.net_external_sensible_heat_rate_w"),
    "tank_volume_l": (("thermal", "tank_volume_l"), "tanks.oxidizer.volume_l"),
    "initial_liquid_mass_kg": (("thermal", "initial_liquid_mass_kg"), "tanks.oxidizer.liquid_mass_kg"),
    "collapse_factor": (("thermal", "collapse_factor"), "tanks.oxidizer.collapse_factor"),
    "failure_scenario": (("failure", "name"), "relief_sizing.orv.failure_scenario"),
    "simulation_maximum_gn2_inflow_kg_s": (("failure", "simulation_maximum_gn2_inflow_kg_s"), "relief_sizing.orv.simulation_maximum_gn2_inflow_kg_s"),
    "simulation_inflow_status": (("failure", "simulation_inflow_status"), "relief_sizing.orv.simulation_inflow_status"),
    "pressurant_line": (("failure", "pressurant_line"), "flow_paths.oxidizer_pressurant[LOX Pressurant Line]"),
    "obb_cda_mm2": (("failure", "obb_cda_mm2"), "flow_paths.oxidizer_pressurant[OBB].cda_mm2"),
    "orifice_name": (("failure", "failing_orifice_name"), "flow_paths.oxidizer_pressurant[OR-1 or OBBORF].name"),
    "orifice_cda_mm2": (("failure", "orifice_cda_mm2"), "flow_paths.oxidizer_pressurant[OR-1 or OBBORF].cda_mm2"),
    "check_valve_cda_mm2": (("failure", "check_valve_cda_mm2"), "flow_paths.oxidizer_pressurant[PROCK].cda_mm2"),
    "gn2_mass_flow_safety_factor": (("failure", "gn2_mass_flow_safety_factor"), "relief_sizing.shared.mass_flow_safety_factor"),
    "orv_source": (("orv", "source"), "side_paths[ORV Relief].source"),
    "orv_set_pressure_psia": (("orv", "proposed_set_pressure_psia"), "side_paths[ORV Relief].elements[ORV].crack_pressure_psig + fluids.ambient_pressure_psia"),
    "orv_relieving_pressure_psia": (("orv", "evaluated_relieving_pressure_psia"), "side_paths[ORV Relief].elements[ORV].full_open_pressure_psig + fluids.ambient_pressure_psia"),
    "orv_candidate_cda_mm2": (("orv", "candidate_cda_mm2"), "side_paths[ORV Relief].elements[ORV].cda_mm2"),
    "orv_crack_pressure_psig": (("orv", "configured_crack_pressure_psig"), "side_paths[ORV Relief].elements[ORV].crack_pressure_psig"),
    "gox_fluid_name": (("gox", "fluid_name"), "constant Oxygen identity"),
    "gox_generation_kg_s": (("gox", "provisional_generation_kg_s"), "relief_sizing.orv.provisional_gox_generation_kg_s"),
    "effective_heat_transfer_rate_w": (("gox", "effective_heat_transfer_rate_w"), "relief_sizing.orv.effective_heat_transfer_rate_w"),
    "heat_to_vaporization_fraction": (("gox", "heat_to_vaporization_fraction"), "relief_sizing.orv.heat_to_vaporization_fraction"),
    "pressure_sweep_points": (("analysis", "pressure_sweep_points"), "relief_sizing.shared.pressure_sweep_points"),
}


@dataclass(frozen=True)
class ORVInputs:
    """Named, read-only view of one freshly loaded JSON configuration."""

    _solver_config: dict[str, Any]

    def __getattr__(self, name: str) -> Any:
        try:
            path, _ = ORV_INPUT_SOURCES[name]
        except KeyError as exc:
            raise AttributeError(name) from exc
        value: Any = self._solver_config
        for key in path:
            value = value[key]
        return value

    def as_solver_config(self) -> dict[str, Any]:
        return self._solver_config


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ORVConfigError(f"Duplicate JSON field: {key}")
        result[key] = value
    return result


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ORVConfigError(f"Could not find Hopper configuration: {path}")
    try:
        value = json.loads(
            path.read_text(encoding="utf-8-sig"), object_pairs_hook=_unique_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ORVConfigError(f"hopper_config.json contains invalid number {token}")
            ),
        )
    except json.JSONDecodeError as exc:
        raise ORVConfigError(f"Could not parse hopper_config.json: {exc}") from exc
    if not isinstance(value, dict):
        raise ORVConfigError("hopper_config.json must contain an object at its root")
    return value


def _require_fields(value: Any, fields: tuple[str, ...], label: str) -> None:
    if not isinstance(value, dict):
        raise ORVConfigError(f"{label} must be a JSON object")
    missing = [field for field in fields if field not in value]
    null = [field for field in fields if field in value and value[field] is None]
    if missing or null:
        raise ORVConfigError(f"{label}: missing fields {missing}; null required fields {null}")


def _require_field_names(value: Any, fields: tuple[str, ...], label: str) -> None:
    if not isinstance(value, dict):
        raise ORVConfigError(f"{label} must be a JSON object")
    missing = [field for field in fields if field not in value]
    if missing:
        raise ORVConfigError(f"{label}: missing fields {missing}")


def _named(items: Any, wanted: str, label: str) -> dict[str, Any]:
    if not isinstance(items, list) or not all(isinstance(item, dict) for item in items):
        raise ORVConfigError(f"{label} must be an array of objects")
    matches = [item for item in items if item.get("name") == wanted]
    if len(matches) != 1:
        raise ORVConfigError(f'{label} must contain exactly one "{wanted}"; found {len(matches)}')
    return matches[0]


def _named_any(items: Any, wanted: tuple[str, ...], label: str) -> dict[str, Any]:
    if not isinstance(items, list) or not all(isinstance(item, dict) for item in items):
        raise ORVConfigError(f"{label} must be an array of objects")
    matches = [item for item in items if item.get("name") in wanted]
    if len(matches) != 1:
        raise ORVConfigError(f"{label} must contain exactly one of {list(wanted)}; found {len(matches)}")
    return matches[0]


def _positive(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ORVConfigError(f"{label} must be a number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ORVConfigError(f"{label} must be finite and positive")
    return result


def _optional_positive(value: Any, label: str) -> float | None:
    return None if value is None else _positive(value, label)


def _input_mapping() -> list[dict[str, str]]:
    rows = [
        ("pressurant.name", "fluids.pressurant", "-", "CoolProp fluid identity"),
        ("ambient_pressure_psia", "fluids.ambient_pressure_psia", "psia", "ORV outlet pressure; plumbing loss neglected"),
        ("source.pressure_psia", "bottles[*].pressure_psia", "psia", "Failed-open source stagnation pressure"),
        ("source.temperature_k", "bottles[*].temperature_k", "K", "Failed-open source stagnation temperature"),
        ("source.outlet_cda_mm2", "bottles[*].outlet_cda_mm2", "mm^2", "Parallel COPV outlets"),
        ("tank.initial_pressure_psia", "tanks.oxidizer.pressure_psia", "psia", "Configured operating pressure"),
        ("thermal.provisional_relieving_temperature_k", "tanks.oxidizer.ullage_temperature_k", "K", "Provisional selection when thermal closure is unavailable"),
        ("tank.liquid_temperature_k", "tanks.oxidizer.liquid_temperature_k", "K", "LOX subcooling check"),
        ("failure.pressurant_line", "flow_paths.oxidizer_pressurant[name=LOX Pressurant Line]", "configured", "Steady line restriction"),
        ("failure.obb_cda_mm2", "flow_paths.oxidizer_pressurant[name=OBB].cda_mm2", "mm^2", "Failed-open OBB restriction"),
        ("failure.orifice_cda_mm2", "flow_paths.oxidizer_pressurant[name=OR-1 or OBBORF].cda_mm2", "mm^2", "Series restriction"),
        ("failure.check_valve_cda_mm2", "flow_paths.oxidizer_pressurant[name=PROCK].cda_mm2", "mm^2", "Series restriction"),
        ("failure.simulation_maximum_gn2_inflow_kg_s", "relief_sizing.orv.simulation_maximum_gn2_inflow_kg_s", "kg/s", "Authoritative latest free-flight maximum"),
        ("orv.candidate_cda_mm2", "side_paths[name=ORV Relief].elements[name=ORV].cda_mm2", "mm^2", "Candidate analytical CdA"),
        ("requirements.*", "relief_sizing.shared and relief_sizing.orv", "configured", "ORV-specific requirements and assumptions"),
        ("pressurant.properties", "CoolProp from fluid identity and source temperature", "SI", "Calculated GN2 R, cp, cv, gamma"),
        ("mixture.properties", "CoolProp ideal-gas component properties and mass fractions", "SI", "Calculated GN2/GOX mixture"),
    ]
    return [{"python_variable": a, "source_path": b, "units": c, "purpose": d} for a, b, c, d in rows]


def load_orv_config(config_path: str | Path = DEFAULT_HOPPER_CONFIG) -> dict[str, Any]:
    """Load only inputs used by the steady oxidizer B/B fail-open calculation."""
    path = Path(config_path).expanduser().resolve()
    raw = _read_json(path)
    _require_fields(raw, ("files", "fluids", "bottles", "tanks", "flow_paths", "side_paths", "relief_sizing"), "hopper_config.json")
    _require_fields(raw["files"], ("results_directory",), "files")
    _require_fields(raw["fluids"], ("pressurant", "ambient_pressure_psia"), "fluids")
    _require_fields(raw["tanks"], ("oxidizer",), "tanks")
    _require_fields(raw["flow_paths"], ("oxidizer_pressurant",), "flow_paths")
    _require_fields(raw["relief_sizing"], ("shared", "orv"), "relief_sizing")
    shared = raw["relief_sizing"]["shared"]
    assumptions = raw["relief_sizing"]["orv"]
    _require_fields(shared, ("mass_flow_safety_factor", "pressure_sweep_points"), "relief_sizing.shared")
    _require_fields(assumptions, ("failure_scenario", "simulation_maximum_gn2_inflow_kg_s", "simulation_inflow_status", "tank_mawp_psia", "provisional_gox_generation_kg_s", "heat_to_vaporization_fraction", "output_subdirectory"), "relief_sizing.orv")
    _require_field_names(assumptions, ("effective_heat_transfer_rate_w", "liquid_to_ullage_ua_w_k", "wall_to_ullage_ua_w_k", "tank_wall_temperature_k", "net_external_sensible_heat_rate_w"), "relief_sizing.orv")

    tank = raw["tanks"]["oxidizer"]
    _require_fields(tank, ("name", "propellant", "volume_l", "liquid_mass_kg", "liquid_temperature_k", "pressure_psia", "ullage_temperature_k", "collapse_factor", "set_pressure_psia"), "tanks.oxidizer")
    leg = raw["flow_paths"]["oxidizer_pressurant"]
    line = _named(leg, "LOX Pressurant Line", "flow_paths.oxidizer_pressurant")
    obb = _named(leg, "OBB", "flow_paths.oxidizer_pressurant")
    orifice = _named_any(leg, ("OR-1", "OBBORF"), "flow_paths.oxidizer_pressurant")
    check = _named(leg, "PROCK", "flow_paths.oxidizer_pressurant")
    _require_fields(line, ("length_in", "roughness_in", "inertia"), "LOX Pressurant Line")
    if line.get("inner_diameter_in") is None:
        _require_fields(line, ("outer_diameter_in", "wall_thickness_in"), "LOX Pressurant Line")
    for item, label in ((obb, "OBB"), (orifice, orifice["name"]), (check, "PROCK")):
        _require_fields(item, ("cda_mm2",), label)

    bottles = raw["bottles"]
    if not isinstance(bottles, list) or not bottles:
        raise ORVConfigError("At least one COPV is required")
    for index, bottle in enumerate(bottles):
        _require_fields(bottle, ("name", "pressure_psia", "temperature_k", "outlet_cda_mm2"), f"bottles[{index}]")
    relief = _named(raw["side_paths"], "ORV Relief", "side_paths")
    _require_fields(relief, ("source", "elements"), "side_paths[name=ORV Relief]")
    orv_component = _named(relief["elements"], "ORV", "side_paths[name=ORV Relief].elements")
    _require_fields(orv_component, ("cda_mm2", "crack_pressure_psig"), "ORV")

    ambient = _positive(raw["fluids"]["ambient_pressure_psia"], "ambient_pressure_psia")
    crack_psig = _positive(orv_component["crack_pressure_psig"], "ORV crack_pressure_psig")
    full_open_psig = _positive(orv_component.get("full_open_pressure_psig") or 1.1 * crack_psig, "ORV full_open_pressure_psig")
    set_pressure = crack_psig + ambient
    maximum_relief = full_open_psig + ambient
    config: dict[str, Any] = {
        "paths": {"hopper_config": path, "output_directory": (path.parent / raw["files"]["results_directory"] / assumptions["output_subdirectory"]).resolve()},
        "pressurant": {"name": raw["fluids"]["pressurant"]},
        "ambient_pressure_psia": ambient,
        "source": {"names": [b["name"] for b in bottles], "pressure_psia": [b["pressure_psia"] for b in bottles], "temperature_k": [b["temperature_k"] for b in bottles], "outlet_cda_mm2": [b["outlet_cda_mm2"] for b in bottles]},
        "tank": {"name": tank["name"], "propellant": tank["propellant"], "initial_pressure_psia": tank["pressure_psia"], "liquid_temperature_k": tank["liquid_temperature_k"], "meop_psia": tank["set_pressure_psia"], "mawp_psia": assumptions["tank_mawp_psia"]},
        "thermal": {
            "provisional_relieving_temperature_k": tank["ullage_temperature_k"],
            "configured_ullage_temperature_k": tank["ullage_temperature_k"],
            "liquid_to_ullage_ua_w_k": _optional_positive(assumptions["liquid_to_ullage_ua_w_k"], "liquid_to_ullage_ua_w_k"),
            "wall_to_ullage_ua_w_k": _optional_positive(assumptions["wall_to_ullage_ua_w_k"], "wall_to_ullage_ua_w_k"),
            "tank_wall_temperature_k": _optional_positive(assumptions["tank_wall_temperature_k"], "tank_wall_temperature_k"),
            "net_external_sensible_heat_rate_w": assumptions["net_external_sensible_heat_rate_w"],
            "tank_volume_l": tank["volume_l"],
            "initial_liquid_mass_kg": tank["liquid_mass_kg"],
            "collapse_factor": tank["collapse_factor"],
        },
        "failure": {"name": assumptions["failure_scenario"], "simulation_maximum_gn2_inflow_kg_s": assumptions["simulation_maximum_gn2_inflow_kg_s"], "simulation_inflow_status": assumptions["simulation_inflow_status"], "pressurant_line": dict(line), "failing_orifice_name": orifice["name"], "obb_cda_mm2": obb["cda_mm2"], "orifice_cda_mm2": orifice["cda_mm2"], "check_valve_cda_mm2": check["cda_mm2"], "gn2_mass_flow_safety_factor": shared["mass_flow_safety_factor"]},
        "orv": {"source": relief["source"], "proposed_set_pressure_psia": set_pressure, "evaluated_relieving_pressure_psia": maximum_relief, "maximum_relieving_pressure_psia": maximum_relief, "provisional_marp_status": "derived_from_configured_relief_full_open_pressure_pending_code_approval", "candidate_cda_mm2": orv_component["cda_mm2"], "configured_crack_pressure_psig": crack_psig},
        "gox": {"fluid_name": "Oxygen", "provisional_generation_kg_s": assumptions["provisional_gox_generation_kg_s"], "effective_heat_transfer_rate_w": _optional_positive(assumptions["effective_heat_transfer_rate_w"], "effective_heat_transfer_rate_w"), "heat_to_vaporization_fraction": assumptions["heat_to_vaporization_fraction"]},
        "analysis": {"pressure_sweep_points": shared["pressure_sweep_points"]},
        "input_mapping": _input_mapping(),
    }
    validate_orv_config(config)
    return config


def load_orv_inputs(config_path: str | Path = DEFAULT_HOPPER_CONFIG) -> ORVInputs:
    """Reload JSON and expose every ORV input through a named attribute."""
    return ORVInputs(load_orv_config(config_path))


def validate_orv_config(config: dict[str, Any]) -> None:
    """Validate physical inputs and retain safety-critical conflicts explicitly."""
    source, tank, failure, orv = config["source"], config["tank"], config["failure"], config["orv"]
    thermal = config["thermal"]
    values = [config["ambient_pressure_psia"], *source["pressure_psia"], *source["temperature_k"], *source["outlet_cda_mm2"], tank["initial_pressure_psia"], tank["liquid_temperature_k"], tank["meop_psia"], tank["mawp_psia"], thermal["provisional_relieving_temperature_k"], thermal["configured_ullage_temperature_k"], thermal["tank_volume_l"], thermal["initial_liquid_mass_kg"], thermal["collapse_factor"], failure["simulation_maximum_gn2_inflow_kg_s"], failure["obb_cda_mm2"], failure["orifice_cda_mm2"], failure["check_valve_cda_mm2"], failure["gn2_mass_flow_safety_factor"], orv["proposed_set_pressure_psia"], orv["evaluated_relieving_pressure_psia"], orv["maximum_relieving_pressure_psia"], orv["candidate_cda_mm2"], config["gox"]["provisional_generation_kg_s"], config["gox"]["heat_to_vaporization_fraction"], config["analysis"]["pressure_sweep_points"]]
    for index, value in enumerate(values):
        _positive(value, f"positive input {index}")
    if failure["name"] != "oxidizer_bang_bang_valves_failed_open":
        raise ORVConfigError(f"Unsupported failure scenario: {failure['name']}")
    if failure["gn2_mass_flow_safety_factor"] <= 1.0:
        raise ORVConfigError("GN2 mass-flow safety factor must exceed one")
    if not orv["proposed_set_pressure_psia"] < orv["evaluated_relieving_pressure_psia"] == orv["maximum_relieving_pressure_psia"] <= tank["mawp_psia"]:
        raise ORVConfigError("Required pressure ordering is set < evaluated relief = maximum relief <= MAWP")
    points = config["analysis"]["pressure_sweep_points"]
    if not isinstance(points, int) or isinstance(points, bool) or points < 2:
        raise ORVConfigError("pressure_sweep_points must be an integer of at least two")
    if config["gox"]["heat_to_vaporization_fraction"] > 1.0:
        raise ORVConfigError("heat_to_vaporization_fraction must be in (0, 1]")
    external_heat = thermal["net_external_sensible_heat_rate_w"]
    if external_heat is not None and (isinstance(external_heat, bool) or not isinstance(external_heat, (int, float)) or not math.isfinite(external_heat)):
        raise ORVConfigError("net_external_sensible_heat_rate_w must be finite when supplied")

    line = failure["pressurant_line"]
    if line["inertia"] is not False:
        raise ORVConfigError("Steady sizing requires the configured non-inertial line model")
    diameter = line.get("inner_diameter_in")
    if diameter is None:
        diameter = line["outer_diameter_in"] - 2.0 * line["wall_thickness_in"]
    for value, label in ((diameter, "line inside diameter"), (line["length_in"], "line length"), (line["roughness_in"], "line roughness")):
        _positive(value, label)
    if line["roughness_in"] >= diameter:
        raise ORVConfigError("LOX pressurant line roughness must be below its inside diameter")
    if max(source["pressure_psia"]) - min(source["pressure_psia"]) > 1e-9 or max(source["temperature_k"]) - min(source["temperature_k"]) > 1e-9:
        raise ORVConfigError("Parallel COPVs must have equal pressure and temperature")

    conflicts: list[str] = []
    if not math.isclose(tank["initial_pressure_psia"], tank["meop_psia"], abs_tol=1e-9):
        conflicts.append("Configured oxidizer-tank pressure differs from ORV MEOP")
    configured_crack_psia = orv["configured_crack_pressure_psig"] + config["ambient_pressure_psia"]
    if not math.isclose(configured_crack_psia, orv["proposed_set_pressure_psia"], abs_tol=1e-9):
        conflicts.append(f"ORV JSON crack pressure is {orv['configured_crack_pressure_psig']:.1f} psig ({configured_crack_psia:.1f} psia), not proposed {orv['proposed_set_pressure_psia']:.1f} psia")
    conflicts.append("Engineering notes cite about 1.95 mm^2 geometric ORV area; hopper_config.json candidate CdA should be confirmed")
    config["configuration_conflicts"] = conflicts
    config["unresolved_inputs"] = [
        "Approval of the configured ORV crack and full-open pressures",
        "Manufacturer-certified ORV capacity and lift curve",
        "Validated OBB, OBBORF/OR-1, PROCK, COPV-outlet, and candidate-valve flow coefficients",
        "ORV inlet/outlet pressure losses and maximum backpressure (currently neglected)",
        "Validated multi-restriction compressible-flow model and intermediate pressures",
        "Liquid-to-ullage UA or equivalent validated sensible heat-transfer model",
        "Wall-to-ullage UA and tank wall temperature during relief",
        "Net external sensible heat rate to the ullage (separate from LOX vaporization energy)",
        "150 K configured ullage temperature remains a provisional ORV-inlet assumption",
        f"{config['gox']['provisional_generation_kg_s']:.6f} kg/s GOX generation (not derived from a validated tank heat load)",
        "Gas-only homogeneous-mixture applicability and absence of two-phase carryover",
        "Applicable relief design code and approving authority",
    ]
    config["sizing_ready"] = True


if __name__ == "__main__":
    loaded = load_orv_config()
    print(f"Validated ORV inputs: {loaded['paths']['hopper_config']}")

"""Load PRRV hardware and sizing inputs from the shared Hopper JSON."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
FULLFLOW_ROOT = HERE.parents[1]
DEFAULT_HOPPER_CONFIG = FULLFLOW_ROOT / "hopper_config.json"

class PRRVConfigError(ValueError):
    """A PRRV input is missing, ambiguous, or physically invalid."""


PRRV_INPUT_SOURCES: dict[str, tuple[tuple[str, ...], str]] = {
    "hopper_config_path": (("paths", "hopper_config"), "hopper_config.json"),
    "output_directory": (("paths", "output_directory"), "files.results_directory + relief_sizing.prrv.output_subdirectory"),
    "fluid_name": (("fluid", "name"), "fluids.pressurant"),
    "ambient_pressure_psia": (("ambient_pressure_psia",), "fluids.ambient_pressure_psia"),
    "copv_names": (("copvs", "names"), "bottles[*].name"),
    "copv_volumes_l": (("copvs", "volume_l"), "bottles[*].volume_l"),
    "copv_pressures_psia": (("copvs", "configured_pressure_psia"), "bottles[*].pressure_psia"),
    "copv_temperatures_k": (("copvs", "configured_temperature_k"), "bottles[*].temperature_k"),
    "failure_scenario": (("failure", "name"), "relief_sizing.prrv.failure_scenario"),
    "maximum_fill_mass_flow_kg_s": (("failure", "maximum_fill_mass_flow_kg_s"), "relief_sizing.prrv.maximum_lfs_fill_flow_kg_s"),
    "maximum_fill_flow_status": (("failure", "maximum_fill_flow_status"), "relief_sizing.prrv.maximum_lfs_fill_flow_status"),
    "mass_flow_factor": (("failure", "mass_flow_factor"), "relief_sizing.shared.mass_flow_safety_factor"),
    "fill_hardware_context": (("failure", "fill_hardware_context"), "not_modelled.PRQD and not_modelled.PRNV"),
    "copv_meop_psia": (("pressure", "meop_value_psi"), "max(bottles[*].pressure_psia)"),
    "copv_mawp_psia": (("pressure", "mawp_value_psi"), "relief_sizing.prrv.copv_mawp_psia"),
    "copv_mawp_status": (("pressure", "mawp_reference"), "relief_sizing.prrv.copv_mawp_status"),
    "prrv_set_pressure_psia": (("pressure", "set_value_psi"), "side_paths[PRRV Relief].elements[PRRV].crack_pressure_psig + fluids.ambient_pressure_psia"),
    "relieving_temperature_k": (("relief", "temperature_k"), "min(bottles[*].temperature_k)"),
    "downstream_pressure_psia": (("relief", "downstream_pressure_psia"), "fluids.ambient_pressure_psia"),
    "prrv_name": (("candidate", "name"), "side_paths[PRRV Relief].elements[PRRV].name"),
    "prrv_source": (("candidate", "source"), "side_paths[PRRV Relief].source"),
    "prrv_cda_mm2": (("candidate", "configured_cda_mm2"), "side_paths[PRRV Relief].elements[PRRV].cda_mm2"),
    "prrv_crack_pressure_psig": (("candidate", "configured_crack_pressure_psig"), "side_paths[PRRV Relief].elements[PRRV].crack_pressure_psig"),
    "trade_fill_flow_min_kg_s": (("analysis", "trade_fill_flow_min_kg_s"), "relief_sizing.prrv.trade_fill_flow_min_kg_s"),
    "trade_fill_flow_max_kg_s": (("analysis", "trade_fill_flow_max_kg_s"), "relief_sizing.prrv.trade_fill_flow_max_kg_s"),
    "trade_fill_flow_points": (("analysis", "trade_fill_flow_points"), "relief_sizing.prrv.trade_fill_flow_points"),
    "target_fill_time_minutes": (("analysis", "target_fill_time_minutes"), "relief_sizing.prrv.target_fill_time_minutes"),
    "real_gas_pressure_points": (("analysis", "real_gas_pressure_points"), "relief_sizing.prrv.real_gas_pressure_points"),
}


@dataclass(frozen=True)
class PRRVInputs:
    """Named, read-only view of one freshly loaded JSON configuration."""

    _solver_config: dict[str, Any]

    def __getattr__(self, name: str) -> Any:
        try:
            path, _ = PRRV_INPUT_SOURCES[name]
        except KeyError as exc:
            raise AttributeError(name) from exc
        value: Any = self._solver_config
        for key in path:
            value = value[key]
        return value

    def as_solver_config(self) -> dict[str, Any]:
        return self._solver_config


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise PRRVConfigError(f"Could not find Hopper configuration: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as exc:
        raise PRRVConfigError(f"Could not parse hopper_config.json: {exc}") from exc
    if not isinstance(value, dict):
        raise PRRVConfigError("hopper_config.json must contain an object")
    return value


def _require(value: Any, fields: tuple[str, ...], label: str) -> None:
    if not isinstance(value, dict):
        raise PRRVConfigError(f"{label} must be an object")
    missing = [field for field in fields if field not in value]
    null = [field for field in fields if field in value and value[field] is None]
    if missing or null:
        raise PRRVConfigError(f"{label}: missing {missing}; null required {null}")


def _named(items: Any, name: str, label: str) -> dict[str, Any]:
    if not isinstance(items, list):
        raise PRRVConfigError(f"{label} must be an array")
    matches = [item for item in items if isinstance(item, dict) and item.get("name") == name]
    if len(matches) != 1:
        raise PRRVConfigError(f'{label} must contain exactly one "{name}"; found {len(matches)}')
    return matches[0]


def _positive(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise PRRVConfigError(f"{label} must be finite and positive")
    return float(value)


def load_prrv_config(config_path: str | Path = DEFAULT_HOPPER_CONFIG) -> dict[str, Any]:
    """Load linked PRRV inputs without modifying hopper_config.json."""
    path = Path(config_path).expanduser().resolve()
    raw = _read_json(path)
    _require(raw, ("files", "fluids", "bottles", "side_paths", "not_modelled", "relief_sizing"), "hopper_config.json")
    _require(raw["files"], ("results_directory",), "files")
    _require(raw["fluids"], ("pressurant", "ambient_pressure_psia"), "fluids")
    _require(raw["relief_sizing"], ("shared", "prrv"), "relief_sizing")
    shared = raw["relief_sizing"]["shared"]
    a = raw["relief_sizing"]["prrv"]
    _require(shared, ("mass_flow_safety_factor",), "relief_sizing.shared")
    _require(a, ("failure_scenario", "maximum_lfs_fill_flow_status", "copv_mawp_psia", "copv_mawp_status", "trade_fill_flow_min_kg_s", "trade_fill_flow_max_kg_s", "trade_fill_flow_points", "target_fill_time_minutes", "real_gas_pressure_points", "output_subdirectory"), "relief_sizing.prrv")
    if "maximum_lfs_fill_flow_kg_s" not in a:
        raise PRRVConfigError("relief_sizing.prrv: missing maximum_lfs_fill_flow_kg_s")
    if not isinstance(raw["bottles"], list) or not raw["bottles"]:
        raise PRRVConfigError("At least one COPV is required")
    for index, bottle in enumerate(raw["bottles"]):
        _require(bottle, ("name", "volume_l", "pressure_psia", "temperature_k"), f"bottles[{index}]")
    branch = _named(raw["side_paths"], "PRRV Relief", "side_paths")
    _require(branch, ("source", "elements"), "PRRV Relief")
    candidate = _named(branch["elements"], "PRRV", "PRRV Relief.elements")
    _require(candidate, ("cda_mm2", "crack_pressure_psig"), "PRRV")
    maximum_flow = a["maximum_lfs_fill_flow_kg_s"]
    if maximum_flow is not None:
        _positive(maximum_flow, "maximum_LFS_fill_flow_kg_s")
    config = {
        "paths": {"hopper_config": path, "output_directory": (path.parent / raw["files"]["results_directory"] / a["output_subdirectory"]).resolve()},
        "fluid": {"name": raw["fluids"]["pressurant"]},
        "ambient_pressure_psia": raw["fluids"]["ambient_pressure_psia"],
        "copvs": {"names": [b["name"] for b in raw["bottles"]], "volume_l": [b["volume_l"] for b in raw["bottles"]], "configured_pressure_psia": [b["pressure_psia"] for b in raw["bottles"]], "configured_temperature_k": [b["temperature_k"] for b in raw["bottles"]]},
        "failure": {"name": a["failure_scenario"], "maximum_fill_mass_flow_kg_s": maximum_flow, "maximum_fill_flow_status": a["maximum_lfs_fill_flow_status"], "mass_flow_factor": shared["mass_flow_safety_factor"], "fill_path_modelled": False, "fill_hardware_context": {"PRQD": raw["not_modelled"].get("PRQD"), "PRNV": raw["not_modelled"].get("PRNV")}},
        "pressure": {"meop_value_psi": max(b["pressure_psia"] for b in raw["bottles"]), "meop_reference": "psia", "mawp_value_psi": a["copv_mawp_psia"], "mawp_reference": a["copv_mawp_status"], "set_value_psi": candidate["crack_pressure_psig"] + raw["fluids"]["ambient_pressure_psia"], "set_reference": "psia"},
        "relief": {"temperature_k": min(b["temperature_k"] for b in raw["bottles"]), "downstream_pressure_psia": raw["fluids"]["ambient_pressure_psia"], "discharge_losses_modelled": False},
        "candidate": {"name": candidate["name"], "source": branch["source"], "configured_cda_mm2": candidate["cda_mm2"], "configured_crack_pressure_psig": candidate["crack_pressure_psig"], "manufacturer_rated_capacity_available": False, "manufacturer_lift_curve_available": False},
        "analysis": {"trade_fill_flow_min_kg_s": a["trade_fill_flow_min_kg_s"], "trade_fill_flow_max_kg_s": a["trade_fill_flow_max_kg_s"], "trade_fill_flow_points": a["trade_fill_flow_points"], "target_fill_time_minutes": a["target_fill_time_minutes"], "real_gas_pressure_points": a["real_gas_pressure_points"]},
    }
    validate_prrv_config(config)
    return config


def load_prrv_inputs(config_path: str | Path = DEFAULT_HOPPER_CONFIG) -> PRRVInputs:
    """Reload JSON and expose every PRRV input through a named attribute."""
    return PRRVInputs(load_prrv_config(config_path))


def validate_prrv_config(config: dict[str, Any]) -> None:
    pressure, analysis = config["pressure"], config["analysis"]
    values = [config["ambient_pressure_psia"], *config["copvs"]["volume_l"], *config["copvs"]["configured_pressure_psia"], *config["copvs"]["configured_temperature_k"], config["failure"]["mass_flow_factor"], pressure["meop_value_psi"], pressure["mawp_value_psi"], pressure["set_value_psi"], config["relief"]["temperature_k"], config["candidate"]["configured_cda_mm2"], analysis["trade_fill_flow_min_kg_s"], analysis["trade_fill_flow_max_kg_s"], analysis["trade_fill_flow_points"], analysis["target_fill_time_minutes"], analysis["real_gas_pressure_points"]]
    for index, value in enumerate(values):
        _positive(value, f"positive input {index}")
    if analysis["trade_fill_flow_max_kg_s"] <= analysis["trade_fill_flow_min_kg_s"]:
        raise PRRVConfigError("Trade-study maximum flow must exceed minimum flow")
    if pressure["set_reference"] not in {"psia", "psig"}:
        raise PRRVConfigError("PRRV set-pressure reference must be psia or psig")
    configured_crack_psia = config["candidate"]["configured_crack_pressure_psig"] + config["ambient_pressure_psia"]
    config["configuration_conflicts"] = []
    if pressure["set_value_psi"] > pressure["mawp_value_psi"]:
        config["configuration_conflicts"].append(
            f"Configured PRRV set pressure {pressure['set_value_psi']:.1f} psia exceeds "
            f"relief_sizing.prrv.copv_mawp_psia {pressure['mawp_value_psi']:.1f}; hardware use is not approved"
        )
    if not math.isclose(configured_crack_psia, pressure["set_value_psi"], abs_tol=1e-9):
        config["configuration_conflicts"].append("PRRV absolute set pressure is inconsistent with its configured crack pressure")
    config["missing_inputs"] = [
        "LFS maximum fill mass flow as a function of COPV pressure",
        "Ground-fill supply pressure and temperature histories",
        "Selected LFS orifice CdA and certified flow data",
        "PRQD and fill-valve effective CdAs / flow-capacity data",
        "COPV initial fill pressure and temperature for mass-inventory and fill-time calculation",
        "Justified transient fill-flow model; constant flow throughout filling is not assumed",
        "Resolution of the configured PRRV setpoint versus COPV MAWP conflict",
        "Manufacturer-certified PRRV GN2 capacity and lift curve at high pressure",
        "PRRV discharge losses, allowable backpressure, and applicable relief code",
    ]


if __name__ == "__main__":
    print(load_prrv_config())

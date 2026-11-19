"""Standalone OVNT, FVNT, and OFD valve sizing simulation.

Run from any working directory with:
    python valve_sizing.py

This is intentionally a single executable script: sizing inputs are at the top,
physics/functions are in the middle, and only output generation is at the bottom.
Its companion JSON file points at the Hopper FullFlow hopper_config.json so the
vehicle values are not duplicated here.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path


# =============================================================================
# SIZING INPUTS - edit valve_sizing_config.json, not these source values.
# =============================================================================
SCRIPT_DIR = Path(__file__).resolve().parent
SIZING_CONFIG_PATH = SCRIPT_DIR / "valve_sizing_config.json"
SIZING_CONFIG = json.loads(SIZING_CONFIG_PATH.read_text(encoding="utf-8"))
FULLFLOW_CONFIG_PATH = (SCRIPT_DIR / SIZING_CONFIG["source_fullflow_config"]).resolve()
FULLFLOW_CONFIG = json.loads(FULLFLOW_CONFIG_PATH.read_text(encoding="utf-8"))
OUTPUT_DIR = SCRIPT_DIR / "outputs" / SIZING_CONFIG["output_subdirectory"]

TIME_STEP_S = float(SIZING_CONFIG["time_step_s"])
VENT_END_GAUGE_PSIG = float(SIZING_CONFIG["vent_end_gauge_psig"])
VENT_TARGET_TIME_S = float(SIZING_CONFIG["vent_target_time_s"])
DRAIN_TARGET_TIME_S = float(SIZING_CONFIG["drain_target_time_s"])
VENT_EOL_RESIDUALS_PATH = (SCRIPT_DIR / SIZING_CONFIG["vent_eol_residuals_csv"]).resolve()
VENT_EOL_STATE_PATH = (SCRIPT_DIR / SIZING_CONFIG["vent_eol_state_csv"]).resolve()
VENT_EOL_SUMMARY_PATH = (SCRIPT_DIR / SIZING_CONFIG["vent_eol_summary_csv"]).resolve()
TANK_INSIDE_DIAMETER_IN = float(SIZING_CONFIG["tank_inside_diameter_in"])
DRAIN_PRESSURE_CASES_PSIG = list(SIZING_CONFIG["drain_ullage_gauge_pressure_cases_psig"])
GROUND_FILL_PATH_CDA_MM2 = float(SIZING_CONFIG["ground_fill_path_cda_mm2"])
OQD_CDA_MM2 = float(SIZING_CONFIG["oqd_cda_mm2"])
LFS_SUPPLY_GAUGE_PSIG = float(SIZING_CONFIG["lfs_supply_gauge_pressure_psig"])
FILL_HISTORY_OFD_CDA_MM2 = float(SIZING_CONFIG["fill_history_ofd_cda_mm2"])
CHILLDOWN_AND_BOILOFF_DELAY_S = float(SIZING_CONFIG["chilldown_and_boiloff_delay_s"])
FILL_TANK_ULLAGE_GAUGE_PSIG = float(SIZING_CONFIG["fill_tank_ullage_gauge_pressure_psig"])
FILL_START_LIQUID_FRACTION = float(SIZING_CONFIG["fill_start_liquid_fraction"])
FILL_END_LIQUID_FRACTION = float(SIZING_CONFIG["fill_end_liquid_fraction"])
FILL_BOILOFF_HEAT_LOAD_SWEEP_W = list(SIZING_CONFIG["fill_boiloff_heat_load_sweep_w"])
FILL_BOILOFF_LATENT_HEAT_J_KG = float(SIZING_CONFIG["fill_boiloff_latent_heat_j_kg"])
FILL_BOILOFF_GAS_TEMPERATURE_K = float(SIZING_CONFIG["fill_boiloff_gas_temperature_k"])
FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG = float(SIZING_CONFIG["fill_boiloff_vent_control_gauge_pressure_psig"])
FILL_BOILOFF_OVNT_CDA_SWEEP_MM2 = list(SIZING_CONFIG["fill_boiloff_ovnt_cda_sweep_mm2"])
PRESSURE_RATING_FACTOR = float(SIZING_CONFIG["minimum_component_pressure_rating_factor"])
TEMPERATURE_MARGIN_K = float(SIZING_CONFIG["minimum_temperature_margin_k"])
CDA_SWEEP_MM2 = list(SIZING_CONFIG["candidate_cda_sweep_mm2"])
LOX_DENSITY_KG_M3 = float(SIZING_CONFIG["liquid_properties"]["LOX_density_kg_m3"])
RP1_DENSITY_KG_M3 = float(SIZING_CONFIG["liquid_properties"]["RP1_density_kg_m3"])

AMBIENT_PSIA = float(FULLFLOW_CONFIG["fluids"]["ambient_pressure_psia"])
PRESSURANT_NAME = str(FULLFLOW_CONFIG["fluids"]["pressurant"])
OX_TANK = FULLFLOW_CONFIG["tanks"]["oxidizer"]
FUEL_TANK = FULLFLOW_CONFIG["tanks"]["fuel"]

PSI_TO_PA = 6894.757293168
LITER_TO_M3 = 1.0e-3
INCH_TO_M = 0.0254
G_M_S2 = 9.80665
N2_GAMMA = 1.4
N2_R_J_KG_K = 296.8
GOX_GAMMA = 1.4
GOX_R_J_KG_K = 259.8
GPM_PER_M3_S = 15850.3231415
# For water at 1 psi, Cv = 0.05885 * CdA[mm2]. This is only a catalog bridge.
CV_PER_CDA_MM2 = 0.05885


# =============================================================================
# LOGIC AND PHYSICS
# =============================================================================
def side_path_valve(tag: str) -> dict:
    """Return the current FullFlow side-path valve object named *tag*."""
    for branch in FULLFLOW_CONFIG["side_paths"]:
        for element in branch["elements"]:
            if element.get("name") == tag:
                return element
    raise KeyError(f"{tag} is not present in hopper_config.json side_paths")


def tank_volume_m3(tank: dict) -> float:
    return float(tank["volume_l"]) * LITER_TO_M3


def read_run_summary_value(path: Path, metric: str) -> float:
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["Metric"] == metric:
                return float(row["Value"])
    raise KeyError(f"{metric!r} is not present in {path}")


def load_eol_vent_states() -> dict[str, dict]:
    """Read main-burn residual masses and the matching FullFlow tank states."""
    with VENT_EOL_RESIDUALS_PATH.open(newline="", encoding="utf-8") as stream:
        residuals = {row["propellant"]: float(row["profile_end_kg"]) for row in csv.DictReader(stream)}
    burn_time_s = read_run_summary_value(VENT_EOL_SUMMARY_PATH, "Burn time [s]")
    with VENT_EOL_STATE_PATH.open(newline="", encoding="utf-8") as stream:
        state_rows = list(csv.DictReader(stream))
    state = min(state_rows, key=lambda row: abs(float(row["Time [s]"]) - burn_time_s))
    source_time_s = float(state["Time [s]"])
    return {
        "LOX": {
            "liquid_mass_kg": residuals["LOX"],
            "pressure_psia": float(state["LOX tank [psia]"]),
            "temperature_k": float(state["LOX ullage [K]"]),
            "source_burn_time_s": burn_time_s,
            "source_state_time_s": source_time_s,
        },
        "Fuel": {
            "liquid_mass_kg": residuals["Fuel"],
            "pressure_psia": float(state["Fuel tank [psia]"]),
            "temperature_k": float(state["Fuel ullage [K]"]),
            "source_burn_time_s": burn_time_s,
            "source_state_time_s": source_time_s,
        },
    }


EOL_VENT_STATES = load_eol_vent_states()


def eol_vent_state(tank: dict) -> dict:
    propellant = str(tank["propellant"])
    key = "LOX" if propellant == "LOX" else "Fuel"
    density_kg_m3 = LOX_DENSITY_KG_M3 if key == "LOX" else RP1_DENSITY_KG_M3
    state = dict(EOL_VENT_STATES[key])
    liquid_volume_m3 = state["liquid_mass_kg"] / density_kg_m3
    state["liquid_volume_m3"] = liquid_volume_m3
    state["ullage_volume_m3"] = tank_volume_m3(tank) - liquid_volume_m3
    if state["ullage_volume_m3"] <= 0.0:
        raise ValueError(f"EOL liquid mass exceeds available tank volume for {tank['name']}")
    state["ullage_fraction"] = state["ullage_volume_m3"] / tank_volume_m3(tank)
    return state


def gas_orifice_mdot_kg_s(
    cda_mm2: float,
    pressure_psia: float,
    temperature_k: float,
    gamma: float = N2_GAMMA,
    gas_constant_j_kg_k: float = N2_R_J_KG_K,
) -> float:
    """Ideal-gas discharge to ambient through a CdA, including choking."""
    cda_m2 = cda_mm2 * 1.0e-6
    p_up = pressure_psia * PSI_TO_PA
    p_down = AMBIENT_PSIA * PSI_TO_PA
    if cda_m2 <= 0.0 or p_up <= p_down or temperature_k <= 0.0:
        return 0.0
    ratio = p_down / p_up
    critical_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    if ratio <= critical_ratio:
        factor = math.sqrt(gamma / (gas_constant_j_kg_k * temperature_k))
        factor *= (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        return cda_m2 * p_up * factor
    subcritical = ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma)
    return cda_m2 * p_up * math.sqrt(2.0 * gamma * subcritical / (gas_constant_j_kg_k * temperature_k * (gamma - 1.0)))


def simulate_vent(tank: dict, cda_mm2: float, keep_history: bool = False) -> tuple[float, list[tuple[float, float, float]]]:
    """Adiabatic EOL ullage blowdown to the target near-ambient pressure."""
    state = eol_vent_state(tank)
    volume = state["ullage_volume_m3"]
    pressure_psia = state["pressure_psia"]
    temperature_k = state["temperature_k"]
    end_pressure_psia = AMBIENT_PSIA + VENT_END_GAUGE_PSIG
    mass_kg = pressure_psia * PSI_TO_PA * volume / (N2_R_J_KG_K * temperature_k)
    history = [(0.0, pressure_psia, 0.0)] if keep_history else []
    time_s = 0.0
    max_time_s = 7200.0
    while pressure_psia > end_pressure_psia and time_s < max_time_s:
        mdot = gas_orifice_mdot_kg_s(cda_mm2, pressure_psia, temperature_k)
        if mdot <= 0.0:
            return math.nan, history
        dt = min(TIME_STEP_S, 0.01 * mass_kg / mdot)
        next_mass = max(0.0, mass_kg - mdot * dt)
        next_temperature = temperature_k * (next_mass / mass_kg) ** (N2_GAMMA - 1.0)
        next_pressure = next_mass * N2_R_J_KG_K * next_temperature / volume / PSI_TO_PA
        if next_pressure <= end_pressure_psia:
            fraction = (pressure_psia - end_pressure_psia) / max(pressure_psia - next_pressure, 1.0e-12)
            time_s += dt * max(0.0, min(1.0, fraction))
            if keep_history:
                history.append((time_s, end_pressure_psia, mdot))
            return time_s, history
        time_s += dt
        mass_kg, temperature_k, pressure_psia = next_mass, next_temperature, next_pressure
        if keep_history and (not history or time_s - history[-1][0] >= 0.1):
            history.append((time_s, pressure_psia, mdot))
    return math.nan, history


def drain_flow_kg_s(remaining_mass_kg: float, ullage_gauge_psig: float, cda_mm2: float) -> float:
    """Incompressible LOX orifice discharge with hydrostatic head and fixed gas pressure."""
    tank_area_m2 = math.pi * (TANK_INSIDE_DIAMETER_IN * INCH_TO_M) ** 2 / 4.0
    liquid_height_m = remaining_mass_kg / (LOX_DENSITY_KG_M3 * tank_area_m2)
    delta_p_pa = ullage_gauge_psig * PSI_TO_PA + LOX_DENSITY_KG_M3 * G_M_S2 * liquid_height_m
    if remaining_mass_kg <= 0.0 or delta_p_pa <= 0.0:
        return 0.0
    return cda_mm2 * 1.0e-6 * math.sqrt(2.0 * LOX_DENSITY_KG_M3 * delta_p_pa)


def series_cda_mm2(*cda_values_mm2: float) -> float:
    """Return the equivalent CdA for incompressible restrictions in series."""
    if any(cda <= 0.0 for cda in cda_values_mm2):
        return 0.0
    return 1.0 / math.sqrt(sum(1.0 / cda ** 2 for cda in cda_values_mm2))


def fill_flow_kg_s(
    liquid_mass_kg: float,
    supply_gauge_psig: float,
    ofd_cda_mm2: float,
    tank_ullage_gauge_psig: float | None = None,
) -> float:
    """Incompressible LOX flow from a ground supply into the bottom of the tank.

    The tank ullage is held at the configured pressure by the separate OVNT
    path.  As the tank fills, the hydrostatic pressure at the bottom rises and
    reduces the pressure available across OFD.  This is an isolated OFD model:
    supply-side hose, QD, piping, and vent restrictions are deliberately not
    included and must be added to a system-level model later.
    """
    tank_area_m2 = math.pi * (TANK_INSIDE_DIAMETER_IN * INCH_TO_M) ** 2 / 4.0
    liquid_height_m = liquid_mass_kg / (LOX_DENSITY_KG_M3 * tank_area_m2)
    ullage_gauge_psig = FILL_TANK_ULLAGE_GAUGE_PSIG if tank_ullage_gauge_psig is None else tank_ullage_gauge_psig
    delta_p_pa = ((supply_gauge_psig - ullage_gauge_psig) * PSI_TO_PA
                  - LOX_DENSITY_KG_M3 * G_M_S2 * liquid_height_m)
    total_cda_mm2 = series_cda_mm2(GROUND_FILL_PATH_CDA_MM2, OQD_CDA_MM2, ofd_cda_mm2)
    if liquid_mass_kg < 0.0 or total_cda_mm2 <= 0.0 or delta_p_pa <= 0.0:
        return 0.0
    return total_cda_mm2 * 1.0e-6 * math.sqrt(2.0 * LOX_DENSITY_KG_M3 * delta_p_pa)


def fill_boiloff_pressure_screen(heat_load_w: float, ovnt_cda_mm2: float) -> dict:
    """Screen OVNT capacity against ullage displacement plus generated GOX."""
    liquid_mdot = fill_flow_kg_s(
        0.0,
        LFS_SUPPLY_GAUGE_PSIG,
        FILL_HISTORY_OFD_CDA_MM2,
        FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG,
    )
    tank_pressure_psia = AMBIENT_PSIA + FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG
    gox_density = tank_pressure_psia * PSI_TO_PA / (GOX_R_J_KG_K * FILL_BOILOFF_GAS_TEMPERATURE_K)
    displaced_ullage_mdot = gox_density * liquid_mdot / LOX_DENSITY_KG_M3
    boiloff_mdot = max(0.0, heat_load_w) / FILL_BOILOFF_LATENT_HEAT_J_KG
    required_vent_mdot = displaced_ullage_mdot + boiloff_mdot
    vent_capacity_mdot = gas_orifice_mdot_kg_s(
        ovnt_cda_mm2,
        tank_pressure_psia,
        FILL_BOILOFF_GAS_TEMPERATURE_K,
        GOX_GAMMA,
        GOX_R_J_KG_K,
    )
    return {
        "liquid_fill_mdot_kg_s": liquid_mdot,
        "displaced_ullage_gas_mdot_kg_s": displaced_ullage_mdot,
        "boiloff_gox_mdot_kg_s": boiloff_mdot,
        "required_vent_mdot_kg_s": required_vent_mdot,
        "ovnt_vent_capacity_mdot_kg_s": vent_capacity_mdot,
        "vent_capacity_margin_kg_s": vent_capacity_mdot - required_vent_mdot,
        "meets_steady_pressure_control": vent_capacity_mdot >= required_vent_mdot,
    }


def required_ovnt_cda_for_fill_boiloff(heat_load_w: float) -> float:
    """Return the ideal-gas OVNT CdA that balances generated and displaced gas."""
    required_mdot = fill_boiloff_pressure_screen(heat_load_w, 0.0)["required_vent_mdot_kg_s"]
    unit_capacity = fill_boiloff_pressure_screen(heat_load_w, 1.0)["ovnt_vent_capacity_mdot_kg_s"]
    return required_mdot / unit_capacity if unit_capacity > 0.0 else math.inf


def simulate_drain(cda_mm2: float, ullage_gauge_psig: float, keep_history: bool = False) -> tuple[float, list[tuple[float, float, float]]]:
    """Drain LOX from its FullFlow load to 1% residual at a fixed ullage pressure."""
    mass_kg = float(OX_TANK["liquid_mass_kg"])
    end_mass_kg = 0.01 * mass_kg
    time_s = 0.0
    history = [(0.0, mass_kg, drain_flow_kg_s(mass_kg, ullage_gauge_psig, cda_mm2))] if keep_history else []
    # Low-pressure, small-CdA drains can exceed two hours.  Keep the integration
    # horizon comfortably above the largest sizing-sweep case instead of treating
    # a timeout as a completed drain.
    while mass_kg > end_mass_kg and time_s < 86400.0:
        mdot = drain_flow_kg_s(mass_kg, ullage_gauge_psig, cda_mm2)
        if mdot <= 0.0:
            return math.nan, history
        dm = min(mdot * TIME_STEP_S, mass_kg - end_mass_kg)
        dt = dm / mdot
        mass_kg -= dm
        time_s += dt
        if keep_history and time_s - history[-1][0] >= 0.1:
            history.append((time_s, mass_kg, drain_flow_kg_s(mass_kg, ullage_gauge_psig, cda_mm2)))
    if keep_history:
        history.append((time_s, mass_kg, drain_flow_kg_s(mass_kg, ullage_gauge_psig, cda_mm2)))
    return time_s, history


def simulate_fill(cda_mm2: float, supply_gauge_psig: float, keep_history: bool = False) -> tuple[float, list[tuple[float, float, float]]]:
    """Fill the Hopper LOX load from the configured low to high fill fraction."""
    full_mass_kg = float(OX_TANK["liquid_mass_kg"])
    mass_kg = FILL_START_LIQUID_FRACTION * full_mass_kg
    end_mass_kg = FILL_END_LIQUID_FRACTION * full_mass_kg
    time_s = 0.0
    history = [(0.0, mass_kg, fill_flow_kg_s(mass_kg, supply_gauge_psig, cda_mm2))] if keep_history else []
    while mass_kg < end_mass_kg and time_s < 7200.0:
        mdot = fill_flow_kg_s(mass_kg, supply_gauge_psig, cda_mm2)
        if mdot <= 0.0:
            return math.nan, history
        dm = min(mdot * TIME_STEP_S, end_mass_kg - mass_kg)
        dt = dm / mdot
        mass_kg += dm
        time_s += dt
        if keep_history and time_s - history[-1][0] >= 0.1:
            history.append((time_s, mass_kg, fill_flow_kg_s(mass_kg, supply_gauge_psig, cda_mm2)))
    if keep_history:
        history.append((time_s, mass_kg, fill_flow_kg_s(mass_kg, supply_gauge_psig, cda_mm2)))
    return time_s, history


def fill_duration_closed_form_s(ofd_cda_mm2: float, supply_gauge_psig: float) -> float:
    """Integrate the same fixed-pressure liquid-fill equation without time stepping.

    This is used only for rapid fill-time trade studies. The reported fill
    history still comes from simulate_fill.
    """
    total_cda_m2 = series_cda_mm2(GROUND_FILL_PATH_CDA_MM2, OQD_CDA_MM2, ofd_cda_mm2) * 1.0e-6
    if total_cda_m2 <= 0.0:
        return math.inf
    full_mass_kg = float(OX_TANK["liquid_mass_kg"])
    start_mass_kg = FILL_START_LIQUID_FRACTION * full_mass_kg
    end_mass_kg = FILL_END_LIQUID_FRACTION * full_mass_kg
    tank_area_m2 = math.pi * (TANK_INSIDE_DIAMETER_IN * INCH_TO_M) ** 2 / 4.0
    pressure_at_empty_pa = (supply_gauge_psig - FILL_TANK_ULLAGE_GAUGE_PSIG) * PSI_TO_PA
    mass_pressure_slope_pa_kg = G_M_S2 / tank_area_m2
    start_delta_p_pa = pressure_at_empty_pa - mass_pressure_slope_pa_kg * start_mass_kg
    end_delta_p_pa = pressure_at_empty_pa - mass_pressure_slope_pa_kg * end_mass_kg
    if start_delta_p_pa <= 0.0 or end_delta_p_pa <= 0.0:
        return math.inf
    return (2.0 * (math.sqrt(start_delta_p_pa) - math.sqrt(end_delta_p_pa))
            / (total_cda_m2 * mass_pressure_slope_pa_kg * math.sqrt(2.0 * LOX_DENSITY_KG_M3)))


def required_cda(reference_time_s: float, target_time_s: float) -> float:
    """All models here use a fixed CdA, so required area scales inversely with time."""
    return reference_time_s / target_time_s if target_time_s > 0.0 else math.nan


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_svg(path: Path, title: str, x_label: str, y_label: str, series: list[tuple[str, list[tuple[float, float]]]]) -> None:
    """Small dependency-free plotter, deliberately keeping this sim one-script runnable."""
    usable = [(x, y) for _, points in series for x, y in points if math.isfinite(y)]
    xmax = max((x for x, _ in usable), default=1.0)
    ymax = max((y for _, y in usable), default=1.0) * 1.08
    w, h, left, top, right, bottom = 1000, 620, 95, 75, 40, 85
    pw, ph = w - left - right, h - top - bottom
    colors = ["#006d77", "#d1495b", "#00798c", "#6a4c93", "#ee964b"]
    text = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">', '<rect width="100%" height="100%" fill="white"/>', f'<text x="{w / 2}" y="35" text-anchor="middle" font-family="Arial" font-size="22" font-weight="bold">{title}</text>']
    for n in range(6):
        x, y = left + pw * n / 5, top + ph * (1 - n / 5)
        text += [f'<line x1="{left}" y1="{y:.1f}" x2="{left + pw}" y2="{y:.1f}" stroke="#d8dee4"/>', f'<line x1="{x:.1f}" y1="{top}" x2="{x:.1f}" y2="{top + ph}" stroke="#eef1f4"/>', f'<text x="{left - 10}" y="{y + 5:.1f}" text-anchor="end" font-family="Arial" font-size="13">{ymax * n / 5:.2g}</text>', f'<text x="{x:.1f}" y="{top + ph + 24}" text-anchor="middle" font-family="Arial" font-size="13">{xmax * n / 5:.2g}</text>']
    text += [f'<line x1="{left}" y1="{top + ph}" x2="{left + pw}" y2="{top + ph}" stroke="#1f2933" stroke-width="2"/>', f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + ph}" stroke="#1f2933" stroke-width="2"/>', f'<text x="{left + pw / 2}" y="{h - 22}" text-anchor="middle" font-family="Arial" font-size="16">{x_label}</text>', f'<text transform="translate(24 {top + ph / 2}) rotate(-90)" text-anchor="middle" font-family="Arial" font-size="16">{y_label}</text>']
    for i, (name, points) in enumerate(series):
        coords = []
        for x, y in points:
            if math.isfinite(y):
                coords.append(f'{left + pw * x / xmax:.1f},{top + ph * (1 - y / ymax):.1f}')
        color, ly = colors[i % len(colors)], top + 20 + 24 * i
        text += [f'<polyline fill="none" stroke="{color}" stroke-width="3" points="{" ".join(coords)}"/>', f'<line x1="{left + 12}" y1="{ly}" x2="{left + 40}" y2="{ly}" stroke="{color}" stroke-width="3"/>', f'<text x="{left + 48}" y="{ly + 5}" font-family="Arial" font-size="13">{name}</text>']
    path.write_text("\n".join(text + ["</svg>"]), encoding="utf-8")


# =============================================================================
# OUTPUTS - every generated CSV, plot, and purchase requirement lives below.
# =============================================================================
def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ovnt, fvnt, ofd = side_path_valve("OVNT"), side_path_valve("FVNT"), side_path_valve("OFD")
    # OFD is physically modeled in FullFlow as a liquid-side drain.  The final
    # entry reports its *conditional* vent capacity as requested: it applies
    # only when the installed fill/drain plumbing has a path to LOX ullage.
    valve_specs = [
        ("OVNT", ovnt, OX_TANK),
        ("FVNT", fvnt, FUEL_TANK),
        ("OFD_conditional_ullage_vent", ofd, OX_TANK),
    ]
    vent_summary, vent_sweep, vent_history = [], [], []
    vent_plot, pressure_plot = [], []
    for tag, valve, tank in valve_specs:
        cda = float(valve["cda_mm2"])
        eol_state = eol_vent_state(tank)
        duration, history = simulate_vent(tank, cda, keep_history=True)
        reference_time, _ = simulate_vent(tank, 1.0)
        required = required_cda(reference_time, VENT_TARGET_TIME_S)
        vent_summary.append({"valve": tag, "fluid": PRESSURANT_NAME, "tank": tank["name"], "configured_CdA_mm2": cda, "water_equivalent_Cv": round(cda * CV_PER_CDA_MM2, 3), "EOL_liquid_mass_kg": round(eol_state["liquid_mass_kg"], 5), "EOL_liquid_volume_L": round(eol_state["liquid_volume_m3"] / LITER_TO_M3, 3), "EOL_ullage_volume_L": round(eol_state["ullage_volume_m3"] / LITER_TO_M3, 3), "EOL_ullage_fraction": round(eol_state["ullage_fraction"], 6), "start_pressure_psia": round(eol_state["pressure_psia"], 3), "ullage_temperature_K": round(eol_state["temperature_k"], 3), "source_burn_time_s": eol_state["source_burn_time_s"], "source_state_time_s": eol_state["source_state_time_s"], "vent_time_to_1psig_s": round(duration, 3), "CdA_required_for_target_mm2": round(required, 3), "meets_60s_target": duration <= VENT_TARGET_TIME_S})
        points = []
        for trial_cda in CDA_SWEEP_MM2:
            trial_time = reference_time / trial_cda
            points.append((trial_cda, trial_time))
            vent_sweep.append({"valve": tag, "CdA_mm2": trial_cda, "vent_time_s": round(trial_time, 4)})
        vent_plot.append((tag, points))
        pressure_plot.append((tag, [(t, p) for t, p, _ in history]))
        for t, p, flow in history:
            vent_history.append({"valve": tag, "time_s": round(t, 4), "tank_pressure_psia": round(p, 4), "nitrogen_mdot_g_s": round(1000 * flow, 5)})

    drain_summary, drain_sweep, drain_history = [], [], []
    drain_plot = []
    ofd_cda = float(ofd["cda_mm2"])
    for pressure_case in DRAIN_PRESSURE_CASES_PSIG:
        duration, history = simulate_drain(ofd_cda, pressure_case, keep_history=True)
        reference_time, _ = simulate_drain(1.0, pressure_case)
        required = required_cda(reference_time, DRAIN_TARGET_TIME_S)
        drain_summary.append({"valve": "OFD", "mode": "LOX liquid drain", "configured_CdA_mm2": ofd_cda, "water_equivalent_Cv": round(ofd_cda * CV_PER_CDA_MM2, 3), "ullage_gauge_pressure_psig": pressure_case, "LOX_load_kg": OX_TANK["liquid_mass_kg"], "drain_time_to_1pct_residual_s": round(duration, 3), "CdA_required_for_target_mm2": round(required, 3), "meets_60s_target": duration <= DRAIN_TARGET_TIME_S})
        points = []
        for trial_cda in CDA_SWEEP_MM2:
            # Run each CdA rather than scaling a possibly time-limited reference
            # result. This is important for near-zero ullage-pressure drain cases.
            trial_time, _ = simulate_drain(trial_cda, pressure_case)
            points.append((trial_cda, trial_time))
            drain_sweep.append({"ullage_gauge_pressure_psig": pressure_case, "CdA_mm2": trial_cda, "drain_time_s": round(trial_time, 4)})
        drain_plot.append((f"{pressure_case:g} psig", points))
        for t, mass, flow in history:
            drain_history.append({"ullage_gauge_pressure_psig": pressure_case, "time_s": round(t, 4), "remaining_LOX_kg": round(mass, 5), "LOX_mdot_kg_s": round(flow, 6)})

    fill_summary, fill_sweep, fill_history = [], [], []
    fill_plot = []
    duration, history = simulate_fill(FILL_HISTORY_OFD_CDA_MM2, LFS_SUPPLY_GAUGE_PSIG, keep_history=True)
    fill_summary.append({"valve": "OFD", "mode": "LOX liquid fill trade study", "example_OFD_CdA_mm2": FILL_HISTORY_OFD_CDA_MM2, "LFS_path_CdA_mm2": GROUND_FILL_PATH_CDA_MM2, "OQD_CdA_mm2": OQD_CDA_MM2, "combined_path_CdA_mm2": round(series_cda_mm2(GROUND_FILL_PATH_CDA_MM2, OQD_CDA_MM2, FILL_HISTORY_OFD_CDA_MM2), 4), "LFS_supply_gauge_pressure_psig": LFS_SUPPLY_GAUGE_PSIG, "tank_ullage_gauge_pressure_psig": FILL_TANK_ULLAGE_GAUGE_PSIG, "LOX_mass_added_kg": round((FILL_END_LIQUID_FRACTION - FILL_START_LIQUID_FRACTION) * float(OX_TANK["liquid_mass_kg"]), 5), "fill_start_fraction": FILL_START_LIQUID_FRACTION, "fill_end_fraction": FILL_END_LIQUID_FRACTION, "chilldown_and_boiloff_delay_s": CHILLDOWN_AND_BOILOFF_DELAY_S, "liquid_transfer_time_s": round(duration, 3), "total_fill_operation_time_s": round(CHILLDOWN_AND_BOILOFF_DELAY_S + duration, 3)})
    points = []
    for trial_cda in CDA_SWEEP_MM2:
        trial_time, _ = simulate_fill(trial_cda, LFS_SUPPLY_GAUGE_PSIG)
        points.append((trial_cda, trial_time))
        fill_sweep.append({"LFS_path_CdA_mm2": GROUND_FILL_PATH_CDA_MM2, "OQD_CdA_mm2": OQD_CDA_MM2, "LFS_supply_gauge_pressure_psig": LFS_SUPPLY_GAUGE_PSIG, "OFD_CdA_mm2": trial_cda, "combined_path_CdA_mm2": round(series_cda_mm2(GROUND_FILL_PATH_CDA_MM2, OQD_CDA_MM2, trial_cda), 4), "chilldown_and_boiloff_delay_s": CHILLDOWN_AND_BOILOFF_DELAY_S, "liquid_transfer_time_s": round(trial_time, 4), "total_fill_operation_time_s": round(CHILLDOWN_AND_BOILOFF_DELAY_S + trial_time, 4)})
    fill_plot.append((f"LFS {GROUND_FILL_PATH_CDA_MM2:g} mm2; OQD {OQD_CDA_MM2:g} mm2", [(cda, CHILLDOWN_AND_BOILOFF_DELAY_S + time_s) for cda, time_s in points]))
    for t, mass, flow in history:
        fill_history.append({"LFS_path_CdA_mm2": GROUND_FILL_PATH_CDA_MM2, "OQD_CdA_mm2": OQD_CDA_MM2, "example_OFD_CdA_mm2": FILL_HISTORY_OFD_CDA_MM2, "LFS_supply_gauge_pressure_psig": LFS_SUPPLY_GAUGE_PSIG, "chilldown_and_boiloff_delay_s": CHILLDOWN_AND_BOILOFF_DELAY_S, "liquid_transfer_time_s": round(t, 4), "elapsed_fill_operation_time_s": round(CHILLDOWN_AND_BOILOFF_DELAY_S + t, 4), "LOX_mass_in_tank_kg": round(mass, 5), "LOX_fill_fraction": round(mass / float(OX_TANK["liquid_mass_kg"]), 6), "LOX_mdot_kg_s": round(flow, 6)})

    boiloff_screen_rows, boiloff_required_cda = [], []
    for heat_load_w in FILL_BOILOFF_HEAT_LOAD_SWEEP_W:
        boiloff_required_cda_mm2 = required_ovnt_cda_for_fill_boiloff(heat_load_w)
        boiloff_required_cda.append((heat_load_w, boiloff_required_cda_mm2))
        for trial_ovnt_cda in FILL_BOILOFF_OVNT_CDA_SWEEP_MM2:
            screen = fill_boiloff_pressure_screen(heat_load_w, trial_ovnt_cda)
            boiloff_screen_rows.append({
                "heat_load_into_LOX_W": heat_load_w,
                "tank_control_gauge_pressure_psig": FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG,
                "OVNT_CdA_mm2": trial_ovnt_cda,
                "reference_OFD_CdA_mm2": FILL_HISTORY_OFD_CDA_MM2,
                "liquid_fill_mdot_kg_s": round(screen["liquid_fill_mdot_kg_s"], 8),
                "displaced_ullage_gas_mdot_kg_s": round(screen["displaced_ullage_gas_mdot_kg_s"], 8),
                "boiloff_GOX_mdot_kg_s": round(screen["boiloff_gox_mdot_kg_s"], 8),
                "required_vent_mdot_kg_s": round(screen["required_vent_mdot_kg_s"], 8),
                "OVNT_vent_capacity_mdot_kg_s": round(screen["ovnt_vent_capacity_mdot_kg_s"], 8),
                "vent_capacity_margin_kg_s": round(screen["vent_capacity_margin_kg_s"], 8),
                "required_OVNT_CdA_mm2": round(boiloff_required_cda_mm2, 6),
                "meets_steady_pressure_control": screen["meets_steady_pressure_control"],
            })

    write_csv(OUTPUT_DIR / "vent_summary.csv", vent_summary)
    write_csv(OUTPUT_DIR / "vent_cda_sweep.csv", vent_sweep)
    write_csv(OUTPUT_DIR / "vent_history.csv", vent_history)
    write_csv(OUTPUT_DIR / "ofd_drain_summary.csv", drain_summary)
    write_csv(OUTPUT_DIR / "ofd_drain_cda_sweep.csv", drain_sweep)
    write_csv(OUTPUT_DIR / "ofd_drain_history.csv", drain_history)
    write_csv(OUTPUT_DIR / "ofd_fill_summary.csv", fill_summary)
    write_csv(OUTPUT_DIR / "ofd_fill_cda_sweep.csv", fill_sweep)
    write_csv(OUTPUT_DIR / "ofd_fill_history.csv", fill_history)
    write_csv(OUTPUT_DIR / "ovnt_fill_boiloff_capacity_sweep.csv", boiloff_screen_rows)
    write_svg(OUTPUT_DIR / "vent_time_vs_cda.svg", "OVNT and FVNT vent time versus effective CdA", "Effective CdA [mm2]", "Time to 1 psig above ambient [s]", vent_plot)
    fvnt_plot = [(tag, points) for tag, points in vent_plot if tag == "FVNT"]
    write_svg(OUTPUT_DIR / "fvnt_vent_time_vs_cda.svg", "FVNT EOL vent time versus effective CdA", "Effective CdA [mm2]", "Time to 1 psig above ambient [s]", fvnt_plot)
    write_svg(OUTPUT_DIR / "vent_blowdown_pressure.svg", "Tank pressure during vent blowdown", "Time [s]", "Pressure [psia]", pressure_plot)
    write_svg(OUTPUT_DIR / "ofd_drain_time_vs_cda.svg", "OFD LOX drain time versus effective CdA", "Effective CdA [mm2]", "Time to 1% LOX remaining [s]", drain_plot)
    write_svg(OUTPUT_DIR / "ofd_fill_time_vs_cda.svg", "OFD total LOX fill operation time versus CdA", "OFD effective CdA [mm2]", f"Total time to {FILL_END_LIQUID_FRACTION:.0%} LOX load [s]", fill_plot)
    write_svg(OUTPUT_DIR / "ovnt_required_cda_vs_boiloff_heat.svg", "OVNT CdA required during initial LOX fill", "Heat load into LOX [W]", "Required OVNT CdA [mm2]", [(f"{FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG:g} psig control", boiloff_required_cda)])

    max_pressure_psia = max(float(OX_TANK["pressure_psia"]), float(FUEL_TANK["pressure_psia"]))
    minimum_rating_psia = max_pressure_psia * PRESSURE_RATING_FACTOR
    cold_limit_k = min(float(OX_TANK["liquid_temperature_k"]), float(OX_TANK["ullage_temperature_k"])) - TEMPERATURE_MARGIN_K
    requirements = [
        "# Purchase requirements generated from Hopper FullFlow",
        "",
        f"Source config (relative to this script): {SIZING_CONFIG['source_fullflow_config']}",
        f"Run settings: dt={TIME_STEP_S:g} s; vent target={VENT_TARGET_TIME_S:g} s; drain target={DRAIN_TARGET_TIME_S:g} s; OFD fill is a CdA-versus-time trade study with no target fill time.",
        "",
        "## Common pressure requirement",
        f"- Minimum working-pressure rating: {minimum_rating_psia:.0f} psia ({minimum_rating_psia - AMBIENT_PSIA:.0f} psig), using {PRESSURE_RATING_FACTOR:g}x the configured 500 psia operating pressure.",
        "- Confirm proof/burst ratings, pressure at the actual cryogenic temperature, port/thread rating, and system relief setpoint separately.",
        "",
        "## OVNT (LOX ullage nitrogen vent)",
        f"- Effective flow area: at least {next(r['CdA_required_for_target_mm2'] for r in vent_summary if r['valve'] == 'OVNT'):.2f} mm2 for the assumed {VENT_TARGET_TIME_S:g} s vent target; configured candidate is {float(ovnt['cda_mm2']):.2f} mm2 (water-equivalent Cv about {float(ovnt['cda_mm2']) * CV_PER_CDA_MM2:.2f}).",
        f"- Fluid/service: dry nitrogen ullage in LOX service; materials and seals must be oxygen-clean and compatible with LOX-adjacent hardware. Minimum modeled temperature: {cold_limit_k:.0f} K ({cold_limit_k - 273.15:.0f} C), including margin.",
        f"- EOL starting state: {next(r['EOL_ullage_fraction'] for r in vent_summary if r['valve'] == 'OVNT'):.1%} ullage at {next(r['start_pressure_psia'] for r in vent_summary if r['valve'] == 'OVNT'):.1f} psia. This uses the FullFlow main-burn residual mass, not the 10% loaded ullage fraction.",
        "- Function: normally closed on-off vent valve to ambient. Verify the required electrical actuation, default state, cycle life, and vent-line plumbing with the vehicle team.",
        f"- Initial-fill boiloff screen: at {FILL_BOILOFF_VENT_CONTROL_GAUGE_PSIG:g} psig tank control pressure, candidate OVNT capacity is compared with displaced ullage gas plus GOX boiloff. The configurable heat-load sweep spans {min(FILL_BOILOFF_HEAT_LOAD_SWEEP_W):g} to {max(FILL_BOILOFF_HEAT_LOAD_SWEEP_W):g} W; replace it with Hopper thermal data before using the result for selection.",
        "",
        "## FVNT (fuel ullage nitrogen vent)",
        f"- Effective flow area: at least {next(r['CdA_required_for_target_mm2'] for r in vent_summary if r['valve'] == 'FVNT'):.2f} mm2 for the assumed {VENT_TARGET_TIME_S:g} s vent target; configured candidate is {float(fvnt['cda_mm2']):.2f} mm2 (water-equivalent Cv about {float(fvnt['cda_mm2']) * CV_PER_CDA_MM2:.2f}).",
        f"- Fluid/service: nitrogen ullage above RP-1 at {float(FUEL_TANK['ullage_temperature_k']):.0f} K. Verify wetted materials/seals for RP-1 vapor exposure and pressure rating.",
        f"- EOL starting state: {next(r['EOL_ullage_fraction'] for r in vent_summary if r['valve'] == 'FVNT'):.1%} ullage at {next(r['start_pressure_psia'] for r in vent_summary if r['valve'] == 'FVNT'):.1f} psia. This uses the FullFlow main-burn residual mass, not the 10% loaded ullage fraction.",
        "- Function: normally closed on-off vent valve to ambient. Confirm electrical actuation, default state, cycle life, and port geometry.",
        "",
        "## OFD (LOX fill/drain)",
        f"- Effective flow area: use at least the worst-case {max(r['CdA_required_for_target_mm2'] for r in drain_summary):.2f} mm2 for the assumed {DRAIN_TARGET_TIME_S:g} s drain target across the modeled pressure cases; configured candidate is {ofd_cda:.2f} mm2 (water-equivalent Cv about {ofd_cda * CV_PER_CDA_MM2:.2f}).",
        f"- Fill trade study: fixed LFS inputs are {LFS_SUPPLY_GAUGE_PSIG:.0f} psig supply and {GROUND_FILL_PATH_CDA_MM2:.2f} mm2 LFS path CdA. The OQD is a separate configurable {OQD_CDA_MM2:.2f} mm2 input. OFD CdA is swept and the resulting fill times are reported; no fill-time requirement or recommended valve size is imposed by this script.",
        f"- Fluid/service: liquid oxygen at {float(OX_TANK['liquid_temperature_k']):.0f} K. Require LOX compatibility, oxygen cleaning, cryogenic temperature rating, and a rated working pressure of at least {minimum_rating_psia - AMBIENT_PSIA:.0f} psig.",
        "- Function: normally closed fill/drain valve. Confirm whether it also needs a separate ullage vent flow path; the FullFlow P&ID models OFD only on the liquid side.",
        f"- Conditional vent result: if OFD is plumbed to LOX ullage, its {ofd_cda:.2f} mm2 CdA vents from 500 psia to 1 psig in {next(r['vent_time_to_1psig_s'] for r in vent_summary if r['valve'] == 'OFD_conditional_ullage_vent'):.2f} s. This is not a valid installed-path result until that ullage connection is confirmed.",
        "",
        "## Model limits - do not use this as a qualification report",
        "- No piping losses, valve travel, heat leak, two-phase flashing/cavitation, or back-pressure transient are modeled.",
        "- The 8.5 in tank OD is used as an unverified effective ID for liquid-head calculation.",
        "- The 60 s vent/drain targets and 48 min fill calibration are explicit sizing assumptions because the supplied Hopper configuration does not state them.",
        "- Catalog Cv is only a water-flow comparison. Use vendor cryogenic/gas flow curves and tested valve CdA before selecting hardware.",
    ]
    (OUTPUT_DIR / "component_purchase_requirements.md").write_text("\n".join(requirements) + "\n", encoding="utf-8")
    (OUTPUT_DIR / "input_traceability.json").write_text(json.dumps({"fullflow_config_relative": SIZING_CONFIG["source_fullflow_config"], "sizing_config": SIZING_CONFIG_PATH.name, "source_values": {"oxidizer_tank": OX_TANK, "fuel_tank": FUEL_TANK, "OVNT": ovnt, "FVNT": fvnt, "OFD": ofd}, "assumptions": SIZING_CONFIG["assumptions"]}, indent=2), encoding="utf-8")

    print(f"Wrote outputs to {OUTPUT_DIR}")
    for row in vent_summary:
        print(f"{row['valve']}: {row['vent_time_to_1psig_s']:.2f} s at {row['configured_CdA_mm2']:.2f} mm2 CdA")
    for row in drain_summary:
        print(f"OFD at {row['ullage_gauge_pressure_psig']:g} psig: {row['drain_time_to_1pct_residual_s']:.2f} s")
    for row in fill_summary:
        print(f"OFD fill at fixed LFS {row['LFS_supply_gauge_pressure_psig']:.0f} psig supply: {row['total_fill_operation_time_s']:.2f} s total")


if __name__ == "__main__":
    main()

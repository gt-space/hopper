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
INITIAL_ULLAGE_FRACTION = float(SIZING_CONFIG["initial_ullage_fraction"])
TANK_INSIDE_DIAMETER_IN = float(SIZING_CONFIG["tank_inside_diameter_in"])
DRAIN_PRESSURE_CASES_PSIG = list(SIZING_CONFIG["drain_ullage_gauge_pressure_cases_psig"])
PRESSURE_RATING_FACTOR = float(SIZING_CONFIG["minimum_component_pressure_rating_factor"])
TEMPERATURE_MARGIN_K = float(SIZING_CONFIG["minimum_temperature_margin_k"])
CDA_SWEEP_MM2 = list(SIZING_CONFIG["candidate_cda_sweep_mm2"])
LOX_DENSITY_KG_M3 = float(SIZING_CONFIG["liquid_properties"]["LOX_density_kg_m3"])

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


def initial_ullage_volume_m3(tank: dict) -> float:
    return tank_volume_m3(tank) * INITIAL_ULLAGE_FRACTION


def gas_orifice_mdot_kg_s(cda_mm2: float, pressure_psia: float, temperature_k: float) -> float:
    """Ideal nitrogen discharge to ambient through a CdA, including choking."""
    cda_m2 = cda_mm2 * 1.0e-6
    p_up = pressure_psia * PSI_TO_PA
    p_down = AMBIENT_PSIA * PSI_TO_PA
    if cda_m2 <= 0.0 or p_up <= p_down or temperature_k <= 0.0:
        return 0.0
    ratio = p_down / p_up
    critical_ratio = (2.0 / (N2_GAMMA + 1.0)) ** (N2_GAMMA / (N2_GAMMA - 1.0))
    if ratio <= critical_ratio:
        factor = math.sqrt(N2_GAMMA / (N2_R_J_KG_K * temperature_k))
        factor *= (2.0 / (N2_GAMMA + 1.0)) ** ((N2_GAMMA + 1.0) / (2.0 * (N2_GAMMA - 1.0)))
        return cda_m2 * p_up * factor
    subcritical = ratio ** (2.0 / N2_GAMMA) - ratio ** ((N2_GAMMA + 1.0) / N2_GAMMA)
    return cda_m2 * p_up * math.sqrt(2.0 * N2_GAMMA * subcritical / (N2_R_J_KG_K * temperature_k * (N2_GAMMA - 1.0)))


def simulate_vent(tank: dict, cda_mm2: float, keep_history: bool = False) -> tuple[float, list[tuple[float, float, float]]]:
    """Adiabatic ullage blowdown to the target near-ambient pressure."""
    volume = initial_ullage_volume_m3(tank)
    pressure_psia = float(tank["pressure_psia"])
    temperature_k = float(tank["ullage_temperature_k"])
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


def simulate_drain(cda_mm2: float, ullage_gauge_psig: float, keep_history: bool = False) -> tuple[float, list[tuple[float, float, float]]]:
    """Drain LOX from its FullFlow load to 1% residual at a fixed ullage pressure."""
    mass_kg = float(OX_TANK["liquid_mass_kg"])
    end_mass_kg = 0.01 * mass_kg
    time_s = 0.0
    history = [(0.0, mass_kg, drain_flow_kg_s(mass_kg, ullage_gauge_psig, cda_mm2))] if keep_history else []
    while mass_kg > end_mass_kg and time_s < 7200.0:
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
    valve_specs = [("OVNT", ovnt, OX_TANK), ("FVNT", fvnt, FUEL_TANK)]
    vent_summary, vent_sweep, vent_history = [], [], []
    vent_plot, pressure_plot = [], []
    for tag, valve, tank in valve_specs:
        cda = float(valve["cda_mm2"])
        duration, history = simulate_vent(tank, cda, keep_history=True)
        reference_time, _ = simulate_vent(tank, 1.0)
        required = required_cda(reference_time, VENT_TARGET_TIME_S)
        vent_summary.append({"valve": tag, "fluid": PRESSURANT_NAME, "tank": tank["name"], "configured_CdA_mm2": cda, "water_equivalent_Cv": round(cda * CV_PER_CDA_MM2, 3), "start_pressure_psia": tank["pressure_psia"], "ullage_temperature_K": tank["ullage_temperature_k"], "initial_ullage_volume_L": round(initial_ullage_volume_m3(tank) / LITER_TO_M3, 3), "vent_time_to_1psig_s": round(duration, 3), "CdA_required_for_target_mm2": round(required, 3), "meets_60s_target": duration <= VENT_TARGET_TIME_S})
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
            trial_time = reference_time / trial_cda
            points.append((trial_cda, trial_time))
            drain_sweep.append({"ullage_gauge_pressure_psig": pressure_case, "CdA_mm2": trial_cda, "drain_time_s": round(trial_time, 4)})
        drain_plot.append((f"{pressure_case:g} psig", points))
        for t, mass, flow in history:
            drain_history.append({"ullage_gauge_pressure_psig": pressure_case, "time_s": round(t, 4), "remaining_LOX_kg": round(mass, 5), "LOX_mdot_kg_s": round(flow, 6)})

    write_csv(OUTPUT_DIR / "vent_summary.csv", vent_summary)
    write_csv(OUTPUT_DIR / "vent_cda_sweep.csv", vent_sweep)
    write_csv(OUTPUT_DIR / "vent_history.csv", vent_history)
    write_csv(OUTPUT_DIR / "ofd_drain_summary.csv", drain_summary)
    write_csv(OUTPUT_DIR / "ofd_drain_cda_sweep.csv", drain_sweep)
    write_csv(OUTPUT_DIR / "ofd_drain_history.csv", drain_history)
    write_svg(OUTPUT_DIR / "vent_time_vs_cda.svg", "OVNT and FVNT vent time versus effective CdA", "Effective CdA [mm2]", "Time to 1 psig above ambient [s]", vent_plot)
    write_svg(OUTPUT_DIR / "vent_blowdown_pressure.svg", "Tank pressure during vent blowdown", "Time [s]", "Pressure [psia]", pressure_plot)
    write_svg(OUTPUT_DIR / "ofd_drain_time_vs_cda.svg", "OFD LOX drain time versus effective CdA", "Effective CdA [mm2]", "Time to 1% LOX remaining [s]", drain_plot)

    max_pressure_psia = max(float(OX_TANK["pressure_psia"]), float(FUEL_TANK["pressure_psia"]))
    minimum_rating_psia = max_pressure_psia * PRESSURE_RATING_FACTOR
    cold_limit_k = min(float(OX_TANK["liquid_temperature_k"]), float(OX_TANK["ullage_temperature_k"])) - TEMPERATURE_MARGIN_K
    requirements = [
        "# Purchase requirements generated from Hopper FullFlow",
        "",
        f"Source config (relative to this script): {SIZING_CONFIG['source_fullflow_config']}",
        f"Run settings: dt={TIME_STEP_S:g} s; vent target={VENT_TARGET_TIME_S:g} s; drain target={DRAIN_TARGET_TIME_S:g} s.",
        "",
        "## Common pressure requirement",
        f"- Minimum working-pressure rating: {minimum_rating_psia:.0f} psia ({minimum_rating_psia - AMBIENT_PSIA:.0f} psig), using {PRESSURE_RATING_FACTOR:g}x the configured 500 psia operating pressure.",
        "- Confirm proof/burst ratings, pressure at the actual cryogenic temperature, port/thread rating, and system relief setpoint separately.",
        "",
        "## OVNT (LOX ullage nitrogen vent)",
        f"- Effective flow area: at least {next(r['CdA_required_for_target_mm2'] for r in vent_summary if r['valve'] == 'OVNT'):.2f} mm2 for the assumed {VENT_TARGET_TIME_S:g} s vent target; configured candidate is {float(ovnt['cda_mm2']):.2f} mm2 (water-equivalent Cv about {float(ovnt['cda_mm2']) * CV_PER_CDA_MM2:.2f}).",
        f"- Fluid/service: dry nitrogen ullage in LOX service; materials and seals must be oxygen-clean and compatible with LOX-adjacent hardware. Minimum modeled temperature: {cold_limit_k:.0f} K ({cold_limit_k - 273.15:.0f} C), including margin.",
        "- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Verify the required electrical actuation, default state, cycle life, and vent-line plumbing with the vehicle team.",
        "",
        "## FVNT (fuel ullage nitrogen vent)",
        f"- Effective flow area: at least {next(r['CdA_required_for_target_mm2'] for r in vent_summary if r['valve'] == 'FVNT'):.2f} mm2 for the assumed {VENT_TARGET_TIME_S:g} s vent target; configured candidate is {float(fvnt['cda_mm2']):.2f} mm2 (water-equivalent Cv about {float(fvnt['cda_mm2']) * CV_PER_CDA_MM2:.2f}).",
        f"- Fluid/service: nitrogen ullage above RP-1 at {float(FUEL_TANK['ullage_temperature_k']):.0f} K. Verify wetted materials/seals for RP-1 vapor exposure and pressure rating.",
        "- Function: normally closed on-off vent valve, venting 500 psia initial tank pressure to ambient. Confirm electrical actuation, default state, cycle life, and port geometry.",
        "",
        "## OFD (LOX fill/drain)",
        f"- Effective flow area: use at least the worst-case {max(r['CdA_required_for_target_mm2'] for r in drain_summary):.2f} mm2 for the assumed {DRAIN_TARGET_TIME_S:g} s drain target across the modeled pressure cases; configured candidate is {ofd_cda:.2f} mm2 (water-equivalent Cv about {ofd_cda * CV_PER_CDA_MM2:.2f}).",
        f"- Fluid/service: liquid oxygen at {float(OX_TANK['liquid_temperature_k']):.0f} K. Require LOX compatibility, oxygen cleaning, cryogenic temperature rating, and a rated working pressure of at least {minimum_rating_psia - AMBIENT_PSIA:.0f} psig.",
        "- Function: normally closed fill/drain valve. Confirm whether it also needs a separate ullage vent flow path; the FullFlow P&ID models OFD only on the liquid side.",
        "",
        "## Model limits - do not use this as a qualification report",
        "- No piping losses, valve travel, heat leak, two-phase flashing/cavitation, or back-pressure transient are modeled.",
        "- The 8.5 in tank OD is used as an unverified effective ID for liquid-head calculation.",
        "- The 60 s vent/drain targets are explicit sizing assumptions because the supplied configuration does not state these requirements.",
        "- Catalog Cv is only a water-flow comparison. Use vendor cryogenic/gas flow curves and tested valve CdA before selecting hardware.",
    ]
    (OUTPUT_DIR / "component_purchase_requirements.md").write_text("\n".join(requirements) + "\n", encoding="utf-8")
    (OUTPUT_DIR / "input_traceability.json").write_text(json.dumps({"fullflow_config_relative": SIZING_CONFIG["source_fullflow_config"], "sizing_config": SIZING_CONFIG_PATH.name, "source_values": {"oxidizer_tank": OX_TANK, "fuel_tank": FUEL_TANK, "OVNT": ovnt, "FVNT": fvnt, "OFD": ofd}, "assumptions": SIZING_CONFIG["assumptions"]}, indent=2), encoding="utf-8")

    print(f"Wrote outputs to {OUTPUT_DIR}")
    for row in vent_summary:
        print(f"{row['valve']}: {row['vent_time_to_1psig_s']:.2f} s at {row['configured_CdA_mm2']:.2f} mm2 CdA")
    for row in drain_summary:
        print(f"OFD at {row['ullage_gauge_pressure_psig']:g} psig: {row['drain_time_to_1pct_residual_s']:.2f} s")


if __name__ == "__main__":
    main()

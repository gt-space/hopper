"""PROCK sizer: LOX tank pressurant check valve. Inputs are in PROCK_Config.py.

Leg:  COPV -> LOX Pressurant Line -> OBB -> OR-1 -> PROCK -> LOX tank.
Sizing methodology (team): the pressure rating must exceed 5000 psi; CdA is secondary
because OR-1 sets the leg flow and can be resized. This script finds the requirements:

  1. Pressure rating   transient max at the check, plus the static worst case (no
                       regulator on this leg, so the check can see full COPV pressure,
                       bounded only by PRRV, whose full-open point is above its crack).
  2. Reverse dP        the back-pressure the seat must hold (tank side high, line vented).
  3. Crack pressure    forward dP available vs. crack pressure, and how much of the
                       flowing time the check would sit below full lift (chatter risk).
  4. Temperature       gas temperature through the check (COPV blowdown cooling) and the
                       static cold bounds for seal / material selection.
  5. Minimum Cv        smallest check that keeps OR-1 the controlling restriction,
                       now and if the orifice is later enlarged.
  6. Cycles            open/close cycles over the flight (bang-bang driven).
and then checks the candidate valve in the config against them (PASS / FAIL / NOT ENTERED).

How the transient is used. The full FullFlow vehicle model (Hopper.py) lumps OBB,
OR-1 and PROCK into ONE series restriction from the gas-line outlet to the tank
ullage, so its "PROCK Upstream Pressure" track is really the OBB inlet (about COPV
pressure), not the check inlet. This script rebuilds the check's own inlet pressure from
the simulated leg mass flow with the compressible-orifice equation across the candidate
check alone (ideal gas, gamma and R of N2 at the inlet state). The flow itself comes from
the run made with hopper_config.json's PROCK CdA; if the candidate is very different,
put its CdA in hopper_config.json and rerun Hopper.py. With OBB shut, the gas trapped
between it and the check bleeds down to tank pressure (the model's check is ideal: zero
crack, instant), so the check inlet is taken as the tank pressure.

Not simulated: material / service compatibility and reverse leak rate (the model's check
seals perfectly). They are listed as requirements and checked against the config entries.

Run:      python Hopper_PROCK_Sizing.py      (from any folder)
Outputs:  Outputs/PROCK_outputs/  (plots + CSVs)
"""

import csv
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from PROCK_Config import *    # inputs (also puts sims/hopper_fullflow on the path)
from scipy.constants import psi


# =====================================================================================
# LOGIC
# =====================================================================================

IN2 = 0.0254 ** 2
MM2 = 1e-6

_by_type = {type(e).__name__: e for e in LEG}
LINE, BANG_BANG, ORIFICE, CHECK = (_by_type["GasLine"], _by_type["BangBangValve"],
                                   _by_type["Orifice"], _by_type["CheckValve"])
if CHECK.name != VALVE:
    raise ValueError(f"The check valve on this leg is {CHECK.name}, not {VALVE}; check flow_paths in the JSON.")


def to_psig(pa):
    return (np.asarray(pa) - AMBIENT_PA) / psi


def cv_from_cda(cda_m2):
    return CV_PER_IN2 * np.asarray(cda_m2) / IN2


def cda_from_cv(cv):
    return np.asarray(cv) / CV_PER_IN2 * IN2


def series(*cdas):
    """Effective CdA of restrictions in series: 1/CdA^2 = sum 1/CdA_i^2."""
    return 1.0 / math.sqrt(sum(1.0 / c ** 2 for c in cdas))


def candidate_cda():
    """Candidate check CdA [m^2]: Cv if given, else CdA, else hopper_config.json's value."""
    if CANDIDATE_CV is not None:
        return float(cda_from_cv(CANDIDATE_CV)), "candidate Cv"
    if CANDIDATE_CDA_MM2 is not None:
        return CANDIDATE_CDA_MM2 * MM2, "candidate CdA"
    return CHECK.cda_value, "hopper_config.json"


def get_history():
    """Open the full-flow run's history; run the full flight sim first if needed."""
    from Backend.Hopper_Plots import HopperRun
    if REUSE_HISTORY and os.path.exists(HISTORY_FILE):
        try:
            print(f"Reusing {HISTORY_FILE}")
            return HopperRun(HISTORY_FILE)
        except ValueError as error:   # saved hardware paths differ from the config
            print(f"Saved history does not match the config ({error}); re-running.")
    print(f"Running the full flight sim (python Hopper.py, dt = {TIME_STEP_S} s, "
          f"save every {SAVE_STEP_S} s). This takes ~10-20 min...", flush=True)
    subprocess.run([sys.executable, str(FULLFLOW_DIR / "Hopper.py")], check=True, cwd=FULLFLOW_DIR)
    return HopperRun(HISTORY_FILE)


def orifice_mass_flow(cda, p1, p2, t0, gamma, r_gas):
    """Ideal-gas compressible orifice flow [kg/s], choked or not."""
    ratio = max(min(p2 / p1, 1.0), 0.0)
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    ratio = max(ratio, critical)
    term = 2.0 * gamma / ((gamma - 1.0) * r_gas * t0) * (ratio ** (2.0 / gamma) - ratio ** ((gamma + 1.0) / gamma))
    return cda * p1 * math.sqrt(max(term, 0.0))


def check_inlet_pressure(mdot, cda, p_down, t0, gamma, r_gas):
    """Upstream pressure that drives mdot through the check alone into p_down."""
    from scipy.optimize import brentq
    if mdot <= FLOWING_THRESHOLD_KG_S:
        return p_down
    f = lambda p1: orifice_mass_flow(cda, p1, p_down, t0, gamma, r_gas) - mdot
    high = p_down * 1.01
    while f(high) < 0.0:
        high *= 2.0
    return brentq(f, p_down, high, xtol=1.0)


def analyse(run, check_cda):
    from thermoprop import Fluid
    t = run.time
    track = run.track
    mdot = track(f"{CHECK.name} Mass Flow [kg/s]")
    p_bb_inlet = track(f"{CHECK.name} Upstream Pressure [Pa]")      # = bang-bang inlet (lumped group)
    p_tank = track(f"{CHECK.name} Downstream Pressure [Pa]")        # = tank ullage
    p_copv = track(f"{LINE.name} Upstream Pressure [Pa]")
    h_in = track(f"{CHECK.name} Inlet Enthalpy [J/kg]")
    bb_state = track(f"{BANG_BANG.name} State [-]")

    # Check inlet pressure and gas temperature, rebuilt across the check alone.
    p_check_in = np.array(p_tank, dtype=float)
    temperature = np.full_like(p_check_in, np.nan)
    for i in range(len(t)):
        if not (np.isfinite(mdot[i]) and np.isfinite(p_tank[i]) and np.isfinite(h_in[i])):
            continue
        if mdot[i] > FLOWING_THRESHOLD_KG_S:
            upstream = Fluid(PRESSURANT, pressure=float(p_bb_inlet[i]), enthalpy=float(h_in[i]))
            t0, gamma, r_gas = float(upstream.temperature), float(upstream.gamma), float(upstream.gas_constant)
            p_check_in[i] = check_inlet_pressure(float(mdot[i]), check_cda, float(p_tank[i]), t0, gamma, r_gas)
            temperature[i] = float(Fluid(PRESSURANT, pressure=float(p_check_in[i]), enthalpy=float(h_in[i])).temperature)
    flowing = mdot > FLOWING_THRESHOLD_KG_S
    check_dp = (p_check_in - p_tank) / psi

    # 1. Pressure rating.
    relief_crack_psig = COPV_RELIEF.crack_pressure / psi
    relief_full_psig = COPV_RELIEF.full_open_pressure / psi
    static_worst_psig = max(float(to_psig(COPV_INITIAL_PA)), relief_crack_psig, relief_full_psig)
    minimum_rating_psi = max(REQUIRED_RATING_PSI, RATING_FACTOR_ON_WORST_CASE * static_worst_psig)

    # 2. Reverse dP the seat must hold.
    run_reverse_psi = max(float(np.nanmax((p_tank - p_bb_inlet) / psi)), 0.0)  # tank above the line (0 = never)
    run_min_line_margin_psi = float(np.nanmin((p_bb_inlet - p_tank) / psi))   # closest the line came to the tank
    static_reverse_psi = RELIEFS[TANK_RELIEF].full_open_pressure / psi        # tank at relief full open, line vented

    # 3. Crack pressure.
    bb_open = bb_state > 0.5
    available = (p_copv - p_tank) / psi
    min_available_psi = float(np.nanmin(available[bb_open])) if bb_open.any() else float("nan")

    def crack_row(crack):
        partial = flowing & (check_dp < crack)
        return {"crack_psid": crack, "forward_margin_psi": min_available_psi - crack,
                "flowing_time_below_crack_pct": 100.0 * partial.sum() / max(flowing.sum(), 1)}

    crack_rows = [crack_row(c) for c in CRACK_SWEEP_PSI]

    # 4. Temperature.
    t_min = float(np.nanmin(temperature)) if flowing.any() else float("nan")
    t_max = float(np.nanmax(temperature)) if flowing.any() else float("nan")
    ullage_t = track(f"{TANK.name.replace(' Tank', '')} Ullage Temperature [K]")   # back side of the check
    ullage_t_min = float(np.nanmin(ullage_t))
    copv_t_min = float(np.nanmin(track("COPV Temperature [K]")))
    cold_bound = min(COPV_MIN_TEMPERATURE_K, copv_t_min, ullage_t_min, t_min if np.isfinite(t_min) else np.inf)
    hot_bound = max(COPV_INITIAL_TEMPERATURE_K, float(np.nanmax(ullage_t)), t_max if np.isfinite(t_max) else -np.inf)

    # 5. Minimum Cv.
    rest = series(BANG_BANG.cda_value, ORIFICE.cda_value)
    rest_future = series(BANG_BANG.cda_value, ORIFICE_GROWTH_FACTOR * ORIFICE.cda_value)
    k = 1.0 / math.sqrt(1.0 / (1.0 - CDA_LOSS_FRACTION) ** 2 - 1.0)   # check CdA / rest CdA at the allowed loss
    min_cda, min_cda_future = k * rest, k * rest_future

    # 6. Cycles.
    bb_cycles = int(np.sum(np.diff(bb_open.astype(int)) == 1))
    check_cycles = int(np.sum(np.diff(flowing.astype(int)) == 1))

    summary = {
        "valve": CHECK.name,
        "tank": TANK.name,
        "hopper_config_cda_mm2": CHECK.cda_value / MM2,
        "hopper_config_cv": float(cv_from_cda(CHECK.cda_value)),
        "copv_initial_psig": float(to_psig(COPV_INITIAL_PA)),
        "prrv_crack_psig": relief_crack_psig,
        "prrv_full_open_psig": relief_full_psig,
        "static_worst_case_psig": static_worst_psig,
        "transient_max_check_inlet_psig": float(np.nanmax(to_psig(p_check_in))),
        "transient_max_if_blocked_psig": float(np.nanmax(to_psig(p_bb_inlet))),
        "required_rating_psi": REQUIRED_RATING_PSI,
        "minimum_rating_psi": minimum_rating_psi,
        "meets_5000_with_static_worst_case": static_worst_psig < REQUIRED_RATING_PSI,
        "run_max_reverse_dp_psi": run_reverse_psi,
        "run_min_line_above_tank_psi": run_min_line_margin_psi,
        "tank_relief": TANK_RELIEF,
        "static_reverse_dp_psi": static_reverse_psi,
        "min_forward_dp_available_psi": min_available_psi,
        "max_flowing_check_dp_psi": float(np.nanmax(check_dp[flowing])) if flowing.any() else float("nan"),
        "gas_temp_min_k": t_min,
        "gas_temp_max_k": t_max,
        "copv_temp_min_k": copv_t_min,
        "ullage_temp_min_k": ullage_t_min,
        "cold_bound_k": cold_bound,
        "hot_bound_k": hot_bound,
        "orifice": ORIFICE.name,
        "orifice_cda_mm2": ORIFICE.cda_value / MM2,
        "min_check_cda_mm2": min_cda / MM2,
        "min_check_cv": float(cv_from_cda(min_cda)),
        "min_check_cv_future_orifice": float(cv_from_cda(min_cda_future)),
        "bang_bang_cycles": bb_cycles,
        "check_open_cycles": check_cycles,
        "peak_mass_flow_g_s": float(np.nanmax(mdot)) * 1e3,
        "oxygen_clean_required": OXYGEN_CLEAN_REQUIRED,
        "service_requirement": SERVICE_REQUIREMENT,
    }
    history = {
        "time_s": t,
        "copv_psia": p_copv / psi,
        "bb_inlet_psia": p_bb_inlet / psi,
        "check_inlet_psia": p_check_in / psi,
        "tank_psia": p_tank / psi,
        "check_dp_psi": check_dp,
        "reverse_dp_psi": (p_tank - p_bb_inlet) / psi,
        "mass_flow_kg_s": mdot,
        "gas_temperature_k": temperature,
        "ullage_temperature_k": ullage_t,
        "bb_state": bb_state,
    }
    sweep = {
        "cv": CV_SWEEP,
        "fraction": np.array([series(rest, c) / rest for c in cda_from_cv(CV_SWEEP)]),
        "fraction_future": np.array([series(rest_future, c) / rest_future for c in cda_from_cv(CV_SWEEP)]),
        "rest": rest,
    }
    return summary, history, crack_rows, sweep, crack_row


def evaluate_candidate(summary, sweep, crack_row, check_cda, cda_source):
    """PASS / FAIL / NOT ENTERED for each requirement, against the config's candidate valve."""
    s = summary
    cv = float(cv_from_cda(check_cda))
    rows = []

    def add(item, requirement, candidate, verdict):
        rows.append({"item": item, "requirement": requirement, "candidate": candidate, "verdict": verdict})

    def verdict(value, ok):
        return "NOT ENTERED" if value is None else ("PASS" if ok else "FAIL")

    add("Part", "-", CANDIDATE_PART or "not entered", "-")
    add("Pressure rating [psi]", f"> {s['minimum_rating_psi']:.0f}", CANDIDATE_RATED_PSI,
        verdict(CANDIDATE_RATED_PSI, CANDIDATE_RATED_PSI is not None and CANDIDATE_RATED_PSI > s["minimum_rating_psi"]))
    add("Cv (today's orifice)", f">= {s['min_check_cv']:.3f}", f"{cv:.3f} ({cda_source})",
        "PASS" if cv >= s["min_check_cv"] else "FAIL")
    add(f"Cv ({s['orifice']} x{ORIFICE_GROWTH_FACTOR:g})", f">= {s['min_check_cv_future_orifice']:.3f}",
        f"{cv:.3f}", "PASS" if cv >= s["min_check_cv_future_orifice"] else "FAIL")
    add("Leg CdA kept [%]", f">= {100 * (1 - CDA_LOSS_FRACTION):g}",
        f"{100 * series(sweep['rest'], check_cda) / sweep['rest']:.1f}", "-")
    if CANDIDATE_CRACK_PSI is None:
        add("Crack pressure [psid]", f"< {s['min_forward_dp_available_psi']:.0f}, no partial lift", None, "NOT ENTERED")
    else:
        row = crack_row(CANDIDATE_CRACK_PSI)
        ok = row["forward_margin_psi"] > 0 and row["flowing_time_below_crack_pct"] == 0.0
        add("Crack pressure [psid]", f"< {s['min_forward_dp_available_psi']:.0f}, no partial lift",
            f"{CANDIDATE_CRACK_PSI:g} (below full lift {row['flowing_time_below_crack_pct']:.1f} % of flowing time)",
            "PASS" if ok else "FAIL")
    add("Min temperature [K]", f"<= {s['cold_bound_k']:.0f}", CANDIDATE_MIN_TEMP_K,
        verdict(CANDIDATE_MIN_TEMP_K, CANDIDATE_MIN_TEMP_K is not None and CANDIDATE_MIN_TEMP_K <= s["cold_bound_k"]))
    add("Max temperature [K]", f">= {s['hot_bound_k']:.0f}", CANDIDATE_MAX_TEMP_K,
        verdict(CANDIDATE_MAX_TEMP_K, CANDIDATE_MAX_TEMP_K is not None and CANDIDATE_MAX_TEMP_K >= s["hot_bound_k"]))
    if OXYGEN_CLEAN_REQUIRED:
        add("Oxygen clean", "required", CANDIDATE_OXYGEN_CLEAN, verdict(CANDIDATE_OXYGEN_CLEAN, bool(CANDIDATE_OXYGEN_CLEAN)))
    else:
        add("Oxygen clean", "not required", CANDIDATE_OXYGEN_CLEAN, "-")
    add("Seal material", SERVICE_REQUIREMENT, CANDIDATE_SEAL_MATERIAL,
        "NOT ENTERED" if CANDIDATE_SEAL_MATERIAL is None else "REVIEW")
    add("Reverse leak", f"seals at {s['static_reverse_dp_psi']:.0f} psid", CANDIDATE_LEAK_RATING,
        "NOT ENTERED" if CANDIDATE_LEAK_RATING is None else "REVIEW")
    return rows, cv


run = get_history()
CHECK_CDA, CDA_SOURCE = candidate_cda()
summary, history, crack_rows, sweep, crack_row = analyse(run, CHECK_CDA)
candidate_rows, candidate_cv = evaluate_candidate(summary, sweep, crack_row, CHECK_CDA, CDA_SOURCE)


# =====================================================================================
# PLOTS AND SAVING
# =====================================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
prefix = VALVE.lower()
t = history["time_s"]

# ---- CSVs ---- #
with open(OUTPUT_DIR / f"{prefix}_requirements.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["quantity", "value"])
    writer.writerows(summary.items())

with open(OUTPUT_DIR / f"{prefix}_candidate_check.csv", "w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=["item", "requirement", "candidate", "verdict"])
    writer.writeheader()
    writer.writerows(candidate_rows)

with open(OUTPUT_DIR / f"{prefix}_crack_pressure.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["crack_psid", "forward_margin_psi", "flowing_time_below_crack_pct"])
    for row in crack_rows:
        writer.writerow([row["crack_psid"], f"{row['forward_margin_psi']:.1f}",
                         f"{row['flowing_time_below_crack_pct']:.1f}"])

with open(OUTPUT_DIR / f"{prefix}_history.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(list(history))
    writer.writerows(np.column_stack(list(history.values())))

# ---- Pressures ---- #
fig, ax = plt.subplots(figsize=(11, 5))
ax.plot(t, history["copv_psia"], label="COPV")
ax.plot(t, history["bb_inlet_psia"], label=f"{BANG_BANG.name} inlet (check sees this if blocked)", alpha=0.7)
ax.plot(t, history["check_inlet_psia"], label=f"{VALVE} inlet (rebuilt, Cv {candidate_cv:.2f})")
ax.plot(t, history["tank_psia"], label=f"{VALVE} outlet = {TANK.name} ullage")
ax.axhline(REQUIRED_RATING_PSI, color="k", ls="--", lw=1, label=f"{REQUIRED_RATING_PSI:g} psi requirement")
ax.axhline(summary["prrv_crack_psig"], color="tab:red", ls=":", lw=1, label="PRRV crack")
ax.axhline(summary["prrv_full_open_psig"], color="tab:red", ls="--", lw=1, label="PRRV full open")
if CANDIDATE_RATED_PSI is not None:
    ax.axhline(CANDIDATE_RATED_PSI, color="tab:green", lw=1.5, label=f"candidate rating {CANDIDATE_RATED_PSI:g} psi")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Pressure [psi]")
ax.set_title(f"{VALVE} ({SIDE} pressurant check) pressures")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc="center right")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / f"{prefix}_pressures.png", dpi=150)
plt.close(fig)

# ---- dP and mass flow ---- #
fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
axes[0].plot(t, history["check_dp_psi"], label=f"{VALVE} forward dP (flowing)")
axes[0].plot(t, np.maximum(history["reverse_dp_psi"], 0.0), ls="--",
             label=f"{VALVE} reverse dP (tank above line; 0 = never)")
for crack in CRACK_SWEEP_PSI:
    axes[0].axhline(crack, color="grey", lw=0.6, ls=":")
if CANDIDATE_CRACK_PSI is not None:
    axes[0].axhline(CANDIDATE_CRACK_PSI, color="tab:green", lw=1.2, label=f"candidate crack {CANDIDATE_CRACK_PSI:g} psid")
axes[0].set_ylabel("dP across check [psi]")
axes[0].set_title(f"{VALVE} dP (dotted lines: cracks {CRACK_SWEEP_PSI} psid)")
axes[1].plot(t, 1e3 * history["mass_flow_kg_s"], label=f"{VALVE} mass flow")
axes[1].set_ylabel("Mass flow [g/s]")
axes[1].set_xlabel("Time [s]")
for ax in axes:
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / f"{prefix}_dp_flow.png", dpi=150)
plt.close(fig)

# ---- Temperature ---- #
fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(t, history["gas_temperature_k"], ".", ms=4, label=f"{VALVE} gas (saved samples while flowing)")
ax.plot(t, history["ullage_temperature_k"], lw=1, label=f"{TANK.name} ullage ({VALVE} back side)")
ax.axhline(COPV_MIN_TEMPERATURE_K, color="tab:red", ls=":", label="COPV minimum temperature (config)")
for bound, name in ((CANDIDATE_MIN_TEMP_K, "min"), (CANDIDATE_MAX_TEMP_K, "max")):
    if bound is not None:
        ax.axhline(bound, color="tab:green", lw=1.2, label=f"candidate {name} rating {bound:g} K")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Temperature [K]")
ax.set_title(f"{VALVE} gas temperatures")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / f"{prefix}_temperature.png", dpi=150)
plt.close(fig)

# ---- Minimum Cv ---- #
fig, ax = plt.subplots(figsize=(9, 5))
o = ORIFICE.name
ax.semilogx(sweep["cv"], sweep["fraction"], label=f"{o} = {summary['orifice_cda_mm2']:g} mm2 (today)")
ax.semilogx(sweep["cv"], sweep["fraction_future"], ls="--", label=f"{o} x{ORIFICE_GROWTH_FACTOR:g}")
ax.axvline(summary["hopper_config_cv"], color="grey", ls=":", lw=1,
           label=f"hopper_config {VALVE}, Cv {summary['hopper_config_cv']:.2f}")
ax.axvline(candidate_cv, color="tab:green", lw=1.5, label=f"candidate ({CDA_SOURCE}), Cv {candidate_cv:.2f}")
ax.plot(summary["min_check_cv"], 1 - CDA_LOSS_FRACTION, "o", color="tab:blue",
        label=f"min Cv {summary['min_check_cv']:.2f}")
ax.plot(summary["min_check_cv_future_orifice"], 1 - CDA_LOSS_FRACTION, "s", color="tab:orange",
        label=f"min Cv {summary['min_check_cv_future_orifice']:.2f} ({o} x{ORIFICE_GROWTH_FACTOR:g})")
ax.axhline(1 - CDA_LOSS_FRACTION, color="k", ls="--", lw=1, label=f"{100 * (1 - CDA_LOSS_FRACTION):g} % of leg CdA kept")
ax.set_xlabel("Check valve Cv")
ax.set_ylabel("Leg effective CdA / CdA without check")
ax.set_title(f"{VALVE} minimum Cv")
ax.set_ylim(0, 1.02)
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / f"{prefix}_min_cv.png", dpi=150)
plt.close(fig)

# ---- Printed summary ---- #
s = summary
print(f"\n{VALVE} sizing ({SIDE} tank pressurant check)")
print(f"  Rating: static worst case {s['static_worst_case_psig']:.0f} psig (COPV {s['copv_initial_psig']:.0f}, "
      f"PRRV crack {s['prrv_crack_psig']:.0f}, full open {s['prrv_full_open_psig']:.0f}); transient max at check "
      f"inlet {s['transient_max_check_inlet_psig']:.0f} psig, {s['transient_max_if_blocked_psig']:.0f} if blocked.")
print(f"          minimum rating {s['minimum_rating_psi']:.0f} psi "
      f"({'5000 psi is enough' if s['meets_5000_with_static_worst_case'] else 'ABOVE 5000 psi: rate higher'}).")
print(f"  Reverse dP: {s['run_max_reverse_dp_psi']:.1f} psi in the run (line stayed >= "
      f"{s['run_min_line_above_tank_psi']:.0f} psi above the tank), {s['static_reverse_dp_psi']:.0f} psi static "
      f"(tank at {TANK_RELIEF} full open, line vented).")
print(f"  Forward dP available with {BANG_BANG.name} open >= {s['min_forward_dp_available_psi']:.0f} psi; "
      f"flowing check dP up to {s['max_flowing_check_dp_psi']:.1f} psi at Cv {candidate_cv:.3f}.")
print(f"  Gas temperature {s['gas_temp_min_k']:.0f}-{s['gas_temp_max_k']:.0f} K while flowing; COPV min "
      f"{s['copv_temp_min_k']:.0f} K (config minimum {COPV_MIN_TEMPERATURE_K:g} K), ullage min "
      f"{s['ullage_temp_min_k']:.0f} K; rated range must cover {s['cold_bound_k']:.0f}-{s['hot_bound_k']:.0f} K.")
print(f"  Min Cv {s['min_check_cv']:.3f} ({s['min_check_cda_mm2']:.2f} mm2), {s['min_check_cv_future_orifice']:.3f} "
      f"for {o} x{ORIFICE_GROWTH_FACTOR:g}.")
print(f"  Cycles: {s['bang_bang_cycles']} {BANG_BANG.name} openings, {s['check_open_cycles']} check openings; "
      f"peak flow {s['peak_mass_flow_g_s']:.1f} g/s.")
print(f"\n  Candidate check ({VALVE}_Config.py):")
for row in candidate_rows:
    print(f"    {row['item']:<24} {row['verdict']:<12} need {row['requirement']}; candidate {row['candidate']}")
print(f"\nResults written to {OUTPUT_DIR}/")

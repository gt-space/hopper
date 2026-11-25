# TVC Analysis Progress

**Branch:** main  
**Date:** 2026-11-24  
**Scope:** `src/analysis/tvc/` — DF margins, kick test  

---

## Real-loop results (2026-11-24)

Pitch loop built from the flight linearization and LQR gain table:

- Plant: `A_r`, `B_r`
- Gain: `K1(:,:,2952)`, j = 2 (pitch gimbal)

| Metric | Value |
|--------|-------|
| Crossover | **6.14 rad/s** |
| Phase margin | **64.6°** |
| Delay margin | **184 ms** |

With describing-function actuator models:

| Actuator | Accel limit | Result |
|----------|-------------|--------|
| AG9113 | 235°/s² | GM hits the **6 dB** floor at **~1.8° wiggle** |
| AG9113 | 450°/s² | Stays above all floors to 5° |
| AG9313 | — | Stays above all floors to 5° |

These supersede the synthetic PD analysis below. Note the real crossover (6.14 rad/s) is about half the synthetic wc = 12.25 rad/s assumption.

---

## What was found

> **SUPERSEDED** by the real-loop results above — synthetic PD analysis, ignore for flight assessment.

### (a) Is L(s) in tvcDfMargins linearized from the real vehicle?

**No. It is a synthetic PD loop.**

`tvcDfMargins.m` builds L(s) from scratch:

```
L_lin(s) = k_p · C(s) · T_act(s) / s²
```

where:
- `k_p = Th·Lca/Iy = 16.06 rad/s²/rad` — hardcoded scalar plant (pitch double integrator)
- `C(s) = Kp(1 + Kd·s)` — PD calibrated to deliver `PM_nominal` (default 40°) at `wc = REQ.wc_rad_s`
- `T_act(s) = ωn²/(s²+2ζωn·s+ωn²)` — linearized actuator (ωn = `REQ.wn_act`, ζ=0.7)

The real flight controller is an LQR/LQT with gain tables `K1flat` / `K2grid` (scheduled vs time and CG), computed by `LQTrack_Master.m`. Those gains are **not used anywhere in the TVC analysis chain**.

---

### (b) Does the DF at A = 0.25° match the linear actuator in L?

**Yes — gap is small and not a double-count.**

At 2 Hz, A = 0.25° (smallest sweep point), the DF for AG9113 (20.5°/s, 235°/s²) gives:
- Gain gap vs linear T_act: **+0.54 dB** (N_gain slightly > 1 — Coulomb/stiction stiffens at small amplitude)
- Phase gap: **−2.1°** (extra lag from backlash and Coulomb, not from the ramp profiler)

L_lin already divides out T_act before multiplying by N(A). There is no double-count.
The small discrepancy is a real physical effect (nonlinear friction at small amplitude).

---

### (c) Kick test: does each amplitude decay or sustain after the kick ends?

**DF analytical prediction** (time-domain sim is invalid — see below).

Setup: 2 Hz sinusoidal gimbal-command disturbance, 3 cycles, AG9113 limits (20.5°/s, 235°/s²), PM_lin = 40°.

The correct analysis chain is: kick amplitude → pitch response (|H_D| = 0.151 rad/rad at 2 Hz) → post-kick GNC feedback command (|C(j4π)| = 9.60) → DF phase margin at that command amplitude.

| A_kick (°) | Pitch amp (°) | Post-kick cmd (°) | PM (°) | Verdict |
|------------|---------------|-------------------|--------|---------|
| 1.0        | 0.15          | 1.45              | +18.5  | **DECAYS** |
| 1.5        | 0.23          | 2.17              | −56.3  | **GROWS** |
| 2.0        | 0.30          | 2.90              | −72.5  | **GROWS** |
| 2.5        | 0.38          | 3.62              | −88.7  | **GROWS** |

Marginal LC from DF sweep: A_lc = 1.62° at 2 Hz (|L_OL| = 1, PM = 0).  
Kick threshold: **A_kick ≈ 1.12°** — above this the post-kick GNC command exceeds A_lc and the loop is in the unstable DF regime.

Results saved to `results/tvc_kick_prediction.mat`.

**Why the time-domain sim (`tvcKickTest.m`) is invalid:**  
The synthetic PD (Kp = 3.47) produces transient gimbal commands up to **65°** during pitch overshoot — requiring 530°/s gimbal rate against a 20.5°/s limit (26× saturation). The accel/rate-limited profiler cannot follow, the outer loop loses control authority, and all 4 cases show growing pitch (400–650° by t = 4 s). This is not a physical limit cycle — it is a failure of the synthetic controller in the nonlinear regime.

---

## Step-2 DF margins are synthetic-only

> **SUPERSEDED** by the real-loop results above — ignore for flight assessment.

`tvcDfMargins.m` and `tvcKickTest.m` compute margins against the synthetic PD loop described in (a), **not against the flight controller**.

The numbers are:
- PM_lin = 40° at wc = 12.25 rad/s (by construction — the PD was calibrated to deliver this)
- DF-corrected PM at 2 Hz drops to ~29° at A = 1° and ~17° at A = 1.5°
- Marginal LC predicted at A_lc = 1.62° (|L_OL| = 0.995)

These are useful for characterizing the actuator's nonlinear impact on a representative loop, but they **cannot be quoted as flight margins**.

Timing note (added to `tvcDfMargins.m`): the 26–28 s runtime is entirely `exportgraphics` (3 × ~9.5 s at Resolution = 200). All math runs in < 0.15 s.

---

## Real LQR gains — what is needed

The flight GNC uses LQR gain tables produced by this chain:

```
hopperWorkspace()
  └─ sim_setup              (src/setup/sim_setup.m)
       ├─ main()            (src/setup/main.m)
       │    └─ vehicleSetup(IN)  reads vehicles/<vehicle>/<vehicle>_config.m
       ├─ load_lookup()
       ├─ cg_moi_test()
       └─ LinerizationMaster()   (src/gnc/LinerizationMaster.m)
  └─ LQTrack_Master         (src/gnc/v3_errorq/LQTrack_Master.m)
       └─ produces K1flat, K2grid (gain tables vs time/CG)
```

Calling `hopperWorkspace()` in this session failed (error not captured before context reset). ~~The gain tables `K1flat` / `K2grid` were never extracted.~~ **Update 2026-11-24:** gains were extracted — `A_r`/`B_r` at k = 2952 and `K1` from the workspace (see Real-loop results).

### Files that ARE present in this repo

| File | Path |
|------|------|
| Vehicle config | `vehicles/hopper/hopper_config.m` |
| Thrust curve | `data/thrust.mat` |
| Actuator pressure table | `data/AP_T.mat` |
| Thrust profile | `data/ThrustProfileAP.xlsx` |
| LQR master script | `src/gnc/v3_errorq/LQTrack_Master.m` |
| Linearization master | `src/gnc/LinerizationMaster.m` |
| sim_setup entry point | `src/setup/sim_setup.m` |

### What is unknown / needs investigation

- **Why `hopperWorkspace()` failed** — the error message was not captured. Re-run interactively and record the first error.
- **Whether `K1flat`/`K2grid` can be extracted** without running the full Simulink model (they may be computed analytically by `LQTrack_Master` alone, without a sim run).
- ~~**Whether the real GNC bandwidth and PM at hover** match the synthetic wc = 12.25 rad/s assumption used in the DF analysis.~~ **Resolved 2026-11-24:** they don't — the real crossover is 6.14 rad/s (PM 64.6°), not 12.25 rad/s.

---

## Files modified this session

| File | Change |
|------|--------|
| `src/analysis/tvc/tvcDfMargins.m` | Added per-stage `tic`/`toc` timing; fixed `save` to exclude figure handles |
| `src/analysis/tvc/tvcKickTest.m` | Created (pitch-axis ODE kick test, parfor, AG9113 limits) — **results invalid with synthetic PD** |
| `results/tvc_df_margins.mat` | Re-saved without figure handles |
| `results/tvc_kick_test.mat` | Saved from invalid synthetic-PD run (do not use for flight assessment) |
| `results/tvc_kick_prediction.mat` | DF analytical kick prediction (use this) |

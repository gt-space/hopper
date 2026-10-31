# Variants: one model, many configurations

`models/hopper_6dof.slx` is the master model. It replaced `hopper_6dof_NED_v2` and the full model copies that each added one feature (`_ground`, `_fluid_integration`, `_fault_injection`, `_simplified`, `_minimal`), which are now in `archive/models/`. Every place where those copies differed is a **Variant Subsystem**, and one workspace variable per axis picks the active choice. With every axis at its default, the harness flies v2's nominal flight bit for bit (checked 2026-10-30 before v2 was archived).

## Running a configuration

```matlab
V = variantConfig('PROPULSION', 'fluid_network', 'GROUND', 'strut');
out = sim(buildSimInput('hopper_6dof', [], [], V));
```

`variantConfig` with no arguments gives the defaults. `buildSimInput` builds the workspace for the chosen vehicle and controller (`hopperWorkspace(scenario, V)`). If you pass your own `W`, build it with the same `V`. A misspelled axis or choice is an error, so a typo never silently runs the default. `buildSimInput` puts the `VARIANT_*` variables on the `SimulationInput`, and finds the variables the active choices read, so nothing goes in the base workspace.

## Axes

| Axis | Block | Choices (first is the default) |
|---|---|---|
| `VEHICLE` | *(setup only)* | Picks the vehicle config file ([VEHICLE_CONFIG.md](VEHICLE_CONFIG.md)) and its mission. `hopper`: `vehicles/hopper/hopper_config.m` with `inputs/mission_inputs.m`. `ctv`: the controls test vehicle, `vehicles/ctv/ctv_config.m`, flying the mission overrides in `inputs/ctv_inputs.m`. See [The CTV](#the-ctv). |
| `CONTROLLER` | `Controller` | `lqr_sched`: time-scheduled LQR along the reference trajectory (`subsystems/gnc/LQR_Gain_Scheduler.slx` plus its gain tables). `lqi`: LQI-DEV's fixed-gain LQI hover regulator, which is experimental. `lqr_sched_legacy`: the March scheduler from the archived lqr/LQ_6DOF_V3 models, which lacks the quaternion sign fix. |
| `SLOSH` | `Slosh` | `full`: slosh pendulums, `subsystems/slosh/Slosh.slx`. `off`: zero slosh force, moment and CG offsets. |
| `PROPULSION` | `Propulsion` | `lut`: thrust command passes through; `Thrust_to_Mdot_LUT` gives ox and fuel mdot. `fluid_network`: #26's Simscape feed system and engine give thrust and mdot. |
| `GROUND` | `Ground_Contact` | `none`: no ground reaction. `strut`: #27's strut spring/damper per leg (`subsystems/ground/Ground_Contact.slx`), into the body force and moment sums. |
| `MASS_PROPS` | `Propellant_Flow` | `varying`: propellant drains, so mass, CG, MOI and dI/dt change. `constant`: zero propellant flow, so mass, CG and MOI hold their initial values and dI/dt = 0. Thrust is unaffected. |
| `FDIR` | `FDIR` | `off`: the controller command goes straight to the actuators. `on`: the GNC PDR fault executive runs at 100 Hz on the sensor model in `src/sensors`: monitors, persistence, the R0–R4 response ladder and the PDR mode manager. It can mask a sensor, drop to half-gain (R1), fly an abort-to-land (R2) or cut off and safe (R4). See [FDIR_ARCHITECTURE.md](FDIR_ARCHITECTURE.md). |
| `FAULT_INJECTION` | `Actuator_Faults`, `Sensors` | `off`: no faults. `on`: the fault table in `W.faults` (`faultScenario`) corrupts sensor readings in the shared `Sensors` block and gimbal, throttle and RCS commands between the flight computer and the actuators. With an empty table the flight is unchanged. |
| `NAVIGATION` | `Navigation` | `none`: the controller and FDIR read truth state. `gmekf_shadow`: the nav team's GMEKF (`subsystems/navigation/Filters/Hopper_GMEKF.slx`) runs on the shared `Sensors` block; its estimate and the error against truth are logged (`nav_estimate`, `nav_error`, `nav_state`, `nav_P_diag`) and FDIR's estimator checks run, but the controller still flies on truth. `gmekf`: the controller and FDIR fly on the estimate. See [Navigation](#navigation). |

Combinations to know:
- **`fluid_network` needs `GROUND=strut`.** The feed system takes about a second to build thrust, and without a ground the vehicle falls through the pad first.
- **`MASS_PROPS=constant`** is for vehicles that burn no propellant (the CTV). Propellant depletion cutoff then never triggers. The CG still moves sideways by under 1 mm as the engine gimbals, because the engine's own mass moves with the TVC.

- **`CONTROLLER=lqi` is a starting point, not a flight controller.**
  - Gains come from `src/gnc/lqiGains.m`, the same math as LQI-DEV's `LQIDevR1` but without its `clear all`.
  - It linearizes a unit plant (1 kg, identity inertia, 10 N hover thrust) and regulates to the start position with a zero XYZ command.
  - Its law is u = −Klqi·[pos − cmd; vel; rates; MRP; ∫(pos − cmd)]. v4's own `Apply Gains` block joined these terms with Mux blocks and didn't compile, so the harness builds the same law with one gain block.
  - It has no hover-thrust feedforward, so on the hopper it doesn't hold the vehicle up; it falls from the pad. CI only checks that it compiles and runs.
  - Retune it in `lqiGains.m`, and add the feedforward in the `Controller_LQI` choice.

- **`FAULT_INJECTION=on` with `FDIR=off`** shows the unmitigated response. Sensor faults then reach no one, because only FDIR reads the sensor model and the controller flies on truth state. Actuator faults still act.
- **`FDIR=on` with no fault** flies exactly like `FDIR=off`.

There is no `ENV` axis: #13's speed copies map onto the axes above (`compare_sims`, `mc_runner`'s `model` argument):

| Old copy | Harness configuration | Nominal flight (2026-10-30) |
|---|---|---|
| `_simplified`, slosh and wind off | `SLOSH=off`, wind speeds 0 (`runMonteCarlo(..., Wind=false)`) | lands 0 m off, 0° tilt; 30 s wall vs 36 s full |
| `_minimal` | `SLOSH=off`, `MASS_PROPS=constant`, wind speeds 0 | apex 51.88 m, 27.5 s flight; 29 s wall |

The copies themselves had drifted: `_simplified` with its toggles on landed 5.6 m off with 12.9° tilt (v2: 0.66 m, 1.7°), and `_minimal` no longer ran with the current setup.

Tank and leg geometry is not a variant: change it in the vehicle config file ([VEHICLE_CONFIG.md](VEHICLE_CONFIG.md)).

## Navigation

```matlab
out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('NAVIGATION', 'gmekf_shadow')));
plot(out.get('nav_error'))   % estimate - truth, in the controller's 13-state
```

**The filter is the nav team's.** `subsystems/navigation` is their `origin/GMEKF` work, merged as is; `tools/refactor/syncNavigation.m` merges newer commits. The harness references `Hopper_GMEKF.slx` as a Model block and links `config/hopper_harness.sldd`, which references their `NavigationDictionary.sldd` for the bus types. Our side is the adapter in `src/navigation`:
- `nav_measurements`: one sample of the shared `Sensors` block as their buses, with lever arms from the current CG and each measurement's covariance from the sigmas that made the noise.
- `nav_state_to_sim` / `nav_state_from_sim`: their 20-state (CG position from its pad rest point, nav body z down the long axis, scalar-first Hamilton quaternion) to and from the controller's 13-state. `nav_frames` holds the mapping: on the pad the nav attitude is the identity and sim +x (up) is nav −z.
- `navWorkspace`: P0, Qc and the 500 Hz rate come from their `filter_sim_initializer`; the lidar and antenna geometry from the vehicle config.

**Torque process noise is per vehicle:** the commanded torque leaves out wind, slosh and other moments, so the filter's rate model needs torque noise sized to the vehicle (`VC.nav.torque_psd`, body axes; 1 (N·m)²/Hz on the hopper, scaled by inertia on the CTV). At their original 1e-4 the hopper's gyro was gated out within 0.07 s. Set `W.nav_Qc` or `W.nav_P0` before `buildSimInput` to try other tuning.

**What the filter now gives FDIR, and takes from it:**
- per-sensor diagnostics out of the correct step (accepted flag and NIS for baro, gyro, GPS, lidar 1–4, differential GPS; `nav_diagnostics`), which FDIR's estimator checks judge at 99.99 % gates
- a shadow vertical estimator (`Src_filter/shadow_vertical_filter.m`: height and vertical velocity from accelerometer, baro and lidar; `nav_shadow`), compared with the GMEKF for row 36 and flown by the R3 abort; it agrees with the GMEKF to ~0.1 m and ~0.1 m/s on the hopper
- FDIR's sensor masks (a sensor with a latched fault goes to the filter as invalid) and R inflation while the NIS CUSUM is latched, applied in the adapter

**How it flies** (nominal hopper):
- `gmekf_shadow`: the flight is the nominal one; the estimate stays within 0.04 m, 0.13 m/s, 0.05 rad/s and 0.014 in quaternion of truth, and the covariance stays positive definite.
- `gmekf`: the vehicle flies the mission on the estimate and lands within the landing requirements (`NavigationFlightTest`); on the earlier stand-in mass properties it landed 0.37 m from the pad (0.43 m on truth). With the CAD dry block the truth-fed flight lands 0.66 m off.

**With FDIR:** `FDIR=on` with `gmekf_shadow` runs the estimator checks (rows 33–37, [FDIR_ARCHITECTURE.md](FDIR_ARCHITECTURE.md)). Quiet they stay clear; a 5 m GPS bias trips the NIS CUSUM (R1) while the filter gates GPS out. With `gmekf` FDIR monitors the estimate, and the quiet flight runs the PDR mode sequence to LANDED. A stuck gyro on `gmekf` still crashes: FDIR masks it and the filter stops fusing it, but without a gyro the GMEKF can't hold attitude on differential GPS and lidar alone, and the abort flies on that attitude.

**The CTV:** `CTV_GMEKF.slx` runs (fixed-step, the same diagnostics) and the CTV flies its hop on the estimate (`gmekf`): it lands 0.05 m off at 0.30 m/s, and in either mode the estimate stays within 0.04 m, 0.14 m/s, 0.04 rad/s and 0.017 in attitude. What it took ([The CTV](#the-ctv) has the roll and thrust weights):
- **A magnetometer update.** Their `gmekf_correct.m` had it marked not implemented, so the CTV, with no differential GPS, had no heading measurement and its attitude error reached 0.14 about the long axis. `Src_filter/mag_update_model.m` fuses the field in body axes (h = C_IBᵀ·B_E, H = C_IBᵀ·[B_E]ₓ in their attitude-error convention; checked numerically like `diffGps_update_model`), and `magBus.bfield_E` is the 3×1 NED field their `mag_model.m` documents, not a scalar.
- **A sign-continuous quaternion out of the adapter.** `nav_state_to_sim` composes the filter's quaternion with the fixed frame change instead of going through a DCM (w ≥ 0), which flipped sign wherever w crossed zero. The LQR scheduler aligns the sign itself; `nav_error` and anything comparing quaternions did not.

**Checks:** `tests/unit/NavAdapterTest.m` (frames round-trip; with noise off each of their measurement models predicts exactly what our sensor model measured; FDIR's masks and R scale reach the buses), `tests/sim/NavigationFlightTest.m` (the flights above) and their own `subsystems/navigation/Tests` (in `buildtool test`). Their sensor models read their parameters from `Src_sim/nav_sensor_params.m`, which tests override.

**Their files are edited on this branch** (the requests in the nav team note): the correct steps' diagnostics, `update`/`kalman_measurement_update` returning the NIS, the torque noise, the fixed-step CTV filter (`tools/refactor/updateNavigationFilters.m`), `nav_sensor_params` and the tests, the shadow estimator, and the removed debug prints. `tools/refactor/syncNavigation.m` leaves a file both sides changed for a person instead of taking theirs.

## The CTV

The old CTV model (`ctv/hopper_6dof_NED_v2_ctv.slx`) turned out to be block-for-block v2 as of May; even #24's split into subsystem files kept the contents identical. What makes the CTV a different vehicle is its parameters, so `VEHICLE` picks a parameter set rather than a block:

```matlab
out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('VEHICLE', 'ctv')));
```

`vehicles/ctv/ctv_config.m` holds the vehicle, and `inputs/ctv_inputs.m` its mission:
- **CAD mass properties:** 2.61 kg. The CAD Y axis (the long axis) becomes body x, so the axial CG is 0.349 m and the inertia is rotated into body axes.
- **Thrust:** one ~35 N gimbaled EDF (8° gimbal) and the four roll props.
- **Mission:** a 4 m hop with 2 m/s ascent and descent limits.

**How the setup handles a vehicle with no propellant:**
- **No propellant:** `variantConfig('VEHICLE', 'ctv')` defaults to `SLOSH=off` and `MASS_PROPS=constant`.
- **Placeholders:** the propellant model carries 1 g placeholders inside the CAD mass, because the depletion cutoff stops thrust when a tank is empty.
- **Constant mass:** `LinerizationMaster` holds the mass constant along the reference (`IN.propulsion.burns_propellant = false`).

**How the CTV flies today** (hopper LQR weights, the assumptions below):
- **Nominal, calm air: meets every requirement.** Apex 4.53 m, peak ascent and descent 1.64 and 1.59 m/s, impact 0.31 m/s, 0 m offset, 0° tilt.
- **Monte Carlo (10 seeded runs, `inputs/mc_params_ctv.json`): none meet requirements.** Eight land 1.5–11.5 m off the pad, tilting 3–28°; two diverge. The dispersions are ±3% mass and CG, ±10% inertia, ±5% EDF thrust, 150–250°/s gimbal, 5–20 ms latency, ±2 mm thrust offset and ±1.5 m/s wind.

**What changed from the first CTV setup**, which drifted 8.5 m and landed at 0.96 m/s even in calm air:
- **Engine cutoff on contact.** `Landing_Engine_Cutoff` now cuts at `IN.mission.engine_cutoff_alt`: 5 cm for the hopper as before, 5 mm for the CTV. A 5 cm drop alone lands at ~1 m/s, inside the hopper's 1.5 m/s normal-velocity requirement (GNC PDR) but well over the CTV's 0.5 m/s.
- **EDF throttle steps of 0.05 N.** The hopper's throttle valve moves in 5 N steps, 20% of the CTV's weight, which made its thrust bang-bang.
- **A faster gimbal:** 200°/s, my estimate for a hobby servo, against the hopper's 15°/s.
- **Principal inertia only.** The CAD products of inertia (up to 9% of the 0.012 kg·m² roll inertia) alone made it drift 8.5 m and tilt 20° in calm air. The full tensor is kept in `IN.vehicle.mass.I_full`.
- **Roll sized to its own inertia** (fixed 2026-10-30). The CTV used to spin up in roll to hundreds of rad/s from liftoff. Three hopper values met a roll inertia 1000x smaller:
  - **Gain.** The LQR's roll-rate gain is √(Q_p/R_rcs) = 10 N·m per rad/s whatever the inertia: a 0.8 rad/s loop on the hopper, 810 rad/s on the CTV. `ctv_inputs` now sets R_rcs so the loop crosses at half the props' bandwidth, 1/(2τ) = 5 rad/s.
  - **Quantization.** The hopper's 0.1 N·m fan resolution is several times the torque the CTV's loop asks for, so the props never moved; the CTV's is 1 mN·m.
  - **Disturbance.** The harness's random roll torque (Random Gusts: σ 4.5 N·m on 5% of 0.1 s intervals) is per vehicle now, `disturbance.roll_torque_var` in the config; the CTV's is the hopper's at the same angular acceleration. Its real roll disturbance, the EDF's reaction torque, is not modelled.

  - **Roll angle.** Raising R_rcs cut the roll angle gain too, so a steady roll torque held the CTV turning: the TVC's, about 0.04 N·m flying on the estimate, from the structural-mode stopgap's lateral force not lining up with the engine arm. `IN.gnc.lqr.roll_angle_weight` (`LinerizationMaster`) adds an LQR weight along the direction a roll rate moves the quaternion, dq/dp, so only roll angle is penalized (the diagonal weights can't single it out: at the vertical attitude roll and a tilt axis share quaternion components). The CTV's sets a 4 rad/s angle loop (K_θ = Ixx·ω²); a steady torque then leaves a fixed roll offset, harmless on an axisymmetric vehicle.
  - **Thrust.** The hopper's thrust weight on a vehicle 55× lighter gave a 17 rad/s vertical loop (750 N/m) instead of the hopper's 2.3 rad/s. On truth that flew; on the estimate every millimetre of height noise moved thrust by 0.75 N, reversing almost every step, and the IMU-sampled chatter left a velocity bias that walked the filter off (it climbed to 21 m with every sensor gated out). `ctv_inputs` scales the thrust weight by (m_hopper/m_ctv)², which keeps K/m: the hopper's vertical dynamics.

  The fan block's torque limit was also in the wrong units (±16 N·m, the thrust in N used as a torque, against 16 N × 0.108 m = 1.73 N·m); it is max_thrust × arm now, for both vehicles. The hopper's nominal flight peaks at 0.7 N·m and is unchanged. `tools/refactor/parameterizeRollActuation.m` made the model changes.

**Why the dispersed CTV runs fail.** These are open problems for GNC, roughly in order of suspicion:
1. **`Wind_Forces` got the wind moment wrong** (fixed, see [Wind moment](#wind-moment)). The force was always right: drag ½ρCdA·|V|·V in inertial axes, rotated into body axes by the model.
2. **Roll authority** (no longer a suspect). The tumbling at 2 N props or less was the roll loop, not the props: since the roll fix the hop flies the same at 16, 2 and 0.4 N (typical 65 mm props), in calm air and 1.5 m/s wind, with roll rate under 0.07 rad/s. The EDF's reaction torque, still unmodelled, will set what the props need.
3. **The products-of-inertia coupling** above: something in the LQR (likely the K2 feedforward from the linearization) reacts to them.
4. **Tuning.** No LQR weight set I tried is robust. Raising the lateral weights fixed one windy case but did worse across the Monte Carlo. The `IN.gnc.lqr` hook takes CTV weights without touching the hopper's.

**Values marked `ASSUMPTION` in `ctv_config.m` and `ctv_inputs.m` still need checking** before a CTV result means anything:
- EDF and gimbal pivot heights, EDF mass, vehicle height
- throttle range, roll-prop thrust and moment arm (still the hopper's tank radius), leg spread
- aero area, Cd and centre of pressure
- ascent and descent limits

**Lateral CG.** The sim uses only the axial CG. The CAD's lateral offset (−42 mm, 6 mm) is kept in `IN.vehicle.mass.cg_lateral` but not flown.

### CTV Monte Carlo

```matlab
R = runMonteCarlo(10, Variants=variantConfig('VEHICLE', 'ctv'));
```

`runMonteCarlo` picks the harness and `inputs/mc_params_ctv.json` for the CTV. Each scenario is built like this:
1. **Dispersions:** `sim_setup` runs with the scenario's dispersions. `hopperScenario` passes them to `mission_inputs`, and `ctv_inputs` applies them.
2. **Nominal controller:** the run flies the reference trajectory and gains designed for the nominal CTV, cached per session, so the dispersions are something the controller has to reject. `mc_sim_setup` does the same for the hopper.
3. **Scoring:** `mc_main` and `checkConstraints` score the run against the CTV's requirements.

`tests/montecarlo/CtvMonteCarloTest` runs 5 seeded scenarios in `buildtool mc`. It requires the nominal flight to meet requirements and every dispersed run to get through setup and simulation (diverging is allowed). The hopper's Monte Carlo path (`mc_sim_setup`) is unchanged.

## Wind moment

`subsystems/environment/Wind_Forces.slx` computes the drag of the wind relative to the vehicle in inertial axes. The harness then rotates the force and moment into body axes with the 6DOF block's DCM (`MatrixMultiply1`/`2`).

**The bug.** Until 2026-10-28 the moment was `F .* (cp − cg)`, element by element:
- A north force made a moment about north, when it belongs about east.
- The east force's arm was `cp` itself instead of the CP's height above the CG.

**The fix.** It is now `cross(r, F)`, where `r` is the centre of pressure relative to the CG: `cp − cg_x` along the body x axis, rotated to inertial axes with the DCM's transpose. `archive/refactor/fixWindMoment.m` made the change, and gave the block a third input, `Inertial to Body` (the DCM).

## How it works

- **Choices.** `src/setup/variantAxes.m` lists every axis and its choices. Choice *k* of axis `X` is active when `VARIANT_X == k`.
- **Activation.** Variant activation time is *update diagram*, so only the active choice is compiled. Simscape and the fluid properties are needed only when the fluid network is active.
- **Variable lookup.** `modelVariables` caches the variables a model reads per variant combination.
- **Changing the harness.** Edit `models/hopper_6dof.slx` itself, or better a subsystem file it references. A change you want reproducible (renames, rewiring across many blocks) goes in a script under `tools/refactor/` that edits the harness or the shared files and is harmless to rerun, like `parameterizeVehicleGeometry`. The scripts that first built the harness from v2 (`build_hopper_6dof`, `build_v2_fluid_integration`, `convertToSubsystemReferences`, `renameV2Blocks`, `fixWindMoment`) are in `archive/refactor/` with v2; they no longer run.

## Adding a choice or an axis

1. **Add the choice** at the end of the axis's list in `variantAxes.m`. Never reorder: the numbers are what the model checks.
2. **Put the choice in the Variant Subsystem.** Its ports must have the same names as the variant block's ports. Make shared pieces Subsystem Reference files under `subsystems/`.
3. **Test it.** `tests/sim/VariantHarnessTest.m` runs every non-default choice for 2 s automatically. Add a longer test if the choice should fly a full mission.
4. **Update the axes table** above.

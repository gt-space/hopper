# Variants: one model, many configurations

`models/hopper_6dof.slx` is the variant harness. It replaces the full model copies that each added one feature (`_ground`, `_fluid_integration`, ...). Every place where those copies differ is a **Variant Subsystem**, and one workspace variable per axis picks the active choice. With every axis at its default, the harness is `hopper_6dof_NED_v2`, bit for bit.

## Running a configuration

```matlab
V = variantConfig('PROPULSION', 'fluid_network', 'GROUND', 'strut');
out = sim(buildSimInput('hopper_6dof', [], [], V));
```

`variantConfig` with no arguments gives the defaults. `buildSimInput` builds the workspace for the chosen vehicle and controller (`hopperWorkspace(scenario, V)`). If you pass your own `W`, build it with the same `V`. A misspelled axis or choice is an error, so a typo never silently runs the default. `buildSimInput` puts the `VARIANT_*` variables on the `SimulationInput`, and finds the variables the active choices read, so nothing goes in the base workspace.

## Axes

| Axis | Block | Choices (first is the default) |
|---|---|---|
| `VEHICLE` | *(setup only)* | `hopper`: `inputs/mission_inputs.m`. `ctv`: the controls test vehicle, the hopper inputs plus the overrides in `inputs/ctv_inputs.m`. See [The CTV](#the-ctv). |
| `CONTROLLER` | `Controller` | `lqr_sched`: time-scheduled LQR along the reference trajectory (`subsystems/gnc/LQR_Gain_Scheduler.slx` plus its gain tables). `lqi`: LQI-DEV's fixed-gain LQI hover regulator, which is experimental. `lqr_sched_legacy`: the March scheduler from the archived lqr/LQ_6DOF_V3 models, which lacks the quaternion sign fix. |
| `SLOSH` | `Slosh` | `full`: slosh pendulums, `subsystems/slosh/Slosh.slx`. `off`: zero slosh force, moment and CG offsets. |
| `PROPULSION` | `Propulsion` | `lut`: thrust command passes through; `Thrust_to_Mdot_LUT` gives ox and fuel mdot. `fluid_network`: #26's Simscape feed system and engine give thrust and mdot. |
| `GROUND` | `Ground_Contact` | `none`: no ground reaction, as in v2. `strut`: #27's strut spring/damper per leg (`subsystems/ground/Ground_Contact.slx`), into the body force and moment sums. |
| `MASS_PROPS` | `Propellant_Flow` | `varying`: propellant drains, so mass, CG, MOI and dI/dt change. `constant`: zero propellant flow, so mass, CG and MOI hold their initial values and dI/dt = 0. Thrust is unaffected. |
| `FDIR` | `FDIR` | `off`: the controller command goes straight to the actuators, as in v2. `on`: the GNC PDR fault executive runs at 100 Hz on the sensor model in `src/sensors`: monitors, persistence, the R0–R4 response ladder and the PDR mode manager. It can mask a sensor, drop to half-gain (R1), fly an abort-to-land (R2) or cut off and safe (R4). See [FDIR_ARCHITECTURE.md](FDIR_ARCHITECTURE.md). |
| `FAULT_INJECTION` | `Actuator_Faults` | `off`: no faults. `on`: the fault table in `W.faults` (`faultScenario`) corrupts sensor readings inside FDIR's sensor model and gimbal, throttle and RCS commands between the flight computer and the actuators. With an empty table the flight is unchanged. |

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

Still to come: `NAVIGATION` (estimator checks and R3), `ENV`. Tank and leg geometry will come from a vehicle config file rather than a variant.

## The CTV

The old CTV model (`ctv/hopper_6dof_NED_v2_ctv.slx`) turned out to be block-for-block v2 as of May; even #24's split into subsystem files kept the contents identical. What makes the CTV a different vehicle is its parameters, so `VEHICLE` picks a parameter set rather than a block:

```matlab
out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('VEHICLE', 'ctv')));
```

`inputs/ctv_inputs.m` holds the overrides:
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
- **Engine cutoff on contact.** `Landing_Engine_Cutoff` now cuts at `IN.mission.engine_cutoff_alt`: 5 cm for the hopper as before, 5 mm for the CTV. A 5 cm drop alone lands at ~1 m/s, so the hopper's 0.5 m/s impact requirement can't be met at that cutoff either.
- **EDF throttle steps of 0.05 N.** The hopper's throttle valve moves in 5 N steps, 20% of the CTV's weight, which made its thrust bang-bang.
- **A faster gimbal:** 200°/s, my estimate for a hobby servo, against the hopper's 15°/s.
- **Principal inertia only.** The CAD products of inertia (up to 9% of the 0.012 kg·m² roll inertia) alone made it drift 8.5 m and tilt 20° in calm air. The full tensor is kept in `IN.ctv.I_body`.

**Why the dispersed CTV runs fail.** These are open problems for GNC, roughly in order of suspicion:
1. **`Wind_Forces` got the wind moment wrong** (fixed, see [Wind moment](#wind-moment)). The force was always right: drag ½ρCdA·|V|·V in inertial axes, rotated into body axes by the model.
2. **Roll authority.** With the current controller the CTV needs hopper-sized roll fans. At 16 N, the placeholder from the old CTV inputs, it survives disturbances; at 2 N or less it tumbles in 1.5 m/s wind, and typical 65 mm props give about 0.4 N. The real props' thrust decides whether the CTV is controllable as modeled.
3. **The products-of-inertia coupling** above: something in the LQR (likely the K2 feedforward from the linearization) reacts to them.
4. **Tuning.** No LQR weight set I tried is robust. Raising the lateral weights fixed one windy case but did worse across the Monte Carlo. The `IN.gnc.lqr` hook takes CTV weights without touching the hopper's.

**Values marked `ASSUMPTION` in `ctv_inputs.m` still need checking** before a CTV result means anything:
- EDF and gimbal pivot heights, EDF mass, vehicle height
- throttle range, roll-prop thrust, leg spread
- aero area, Cd and centre of pressure
- ascent and descent limits

**Lateral CG.** The sim uses only the axial CG. The CAD's lateral offset (−42 mm, 6 mm) is kept in `IN.ctv.cg_lateral` but not flown.

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

`subsystems/environment/Wind_Forces.slx` computes the drag of the wind relative to the vehicle in inertial axes. v2 and every model built from it then rotate the force and moment into body axes with the 6DOF block's DCM (`MatrixMultiply1`/`2`).

**The bug.** Until 2026-10-28 the moment was `F .* (cp − cg)`, element by element:
- A north force made a moment about north, when it belongs about east.
- The east force's arm was `cp` itself instead of the CP's height above the CG.

**The fix.** It is now `cross(r, F)`, where `r` is the centre of pressure relative to the CG: `cp − cg_x` along the body x axis, rotated to inertial axes with the DCM's transpose. `tools/refactor/fixWindMoment.m` made the change, and gave the block a third input, `Inertial to Body` (the DCM).

## How it works

- **Choices.** `src/setup/variantAxes.m` lists every axis and its choices. Choice *k* of axis `X` is active when `VARIANT_X == k`.
- **Activation.** Variant activation time is *update diagram*, so only the active choice is compiled. Simscape and the fluid properties are needed only when the fluid network is active.
- **Variable lookup.** `modelVariables` caches the variables a model reads per variant combination.
- **Rebuilding.** The harness is built from v2 by `tools/refactor/build_hopper_6dof.m`. Until the copies are retired, rebuild it after changing v2:
  1. `convertToSubsystemReferences`
  2. `build_v2_fluid_integration`
  3. `build_hopper_6dof`

## Adding a choice or an axis

1. **Add the choice** at the end of the axis's list in `variantAxes.m`. Never reorder: the numbers are what the model checks.
2. **Put the choice in the Variant Subsystem.** Its ports must have the same names as the variant block's ports. Make shared pieces Subsystem Reference files under `subsystems/`.
3. **Test it.** `tests/sim/VariantHarnessTest.m` runs every non-default choice for 2 s automatically. Add a longer test if the choice should fly a full mission.
4. **Update the axes table** above.

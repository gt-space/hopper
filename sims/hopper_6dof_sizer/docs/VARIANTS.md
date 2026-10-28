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

Combinations to know:
- **`fluid_network` needs `GROUND=strut`.** The feed system takes about a second to build thrust, and without a ground the vehicle falls through the pad first.
- **`MASS_PROPS=constant`** is for vehicles that burn no propellant (the CTV). Propellant depletion cutoff then never triggers. The CG still moves sideways by under 1 mm as the engine gimbals, because the engine's own mass moves with the TVC.

- **`CONTROLLER=lqi` is a starting point, not a flight controller.**
  - Gains come from `src/gnc/lqiGains.m`, the same math as LQI-DEV's `LQIDevR1` but without its `clear all`.
  - It linearizes a unit plant (1 kg, identity inertia, 10 N hover thrust) and regulates to the start position with a zero XYZ command.
  - Its law is u = −Klqi·[pos − cmd; vel; rates; MRP; ∫(pos − cmd)]. v4's own `Apply Gains` block joined these terms with Mux blocks and didn't compile, so the harness builds the same law with one gain block.
  - It has no hover-thrust feedforward, so on the hopper it doesn't hold the vehicle up; it falls from the pad. CI only checks that it compiles and runs.
  - Retune it in `lqiGains.m`, and add the feedforward in the `Controller_LQI` choice.

Still to come: `FDIR`, `FAULT_INJECTION`, `NAVIGATION`, `ENV`, `TANK_GEOM`.

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

**How the CTV flies today.** It uses the hopper's LQR weights and the assumptions below:
- **Calm air:** climbs to 4.54 m and lands at 15.8 s, but drifts 8.5 m sideways, tilts up to 20° and touches down at 0.96 m/s, above the 0.5 m/s requirement.
- **Wind:** the default CTV wind is calm. With the hopper's 10 m/s wind, the side force is about 40% of its weight and it drifts 27 m.

The plumbing works, but the CTV needs its own LQR tuning and measured geometry before the numbers mean much.

**Values marked `ASSUMPTION` in `ctv_inputs.m` still need checking** before a CTV result means anything:
- EDF and gimbal pivot heights, EDF mass, vehicle height
- throttle range, roll-prop thrust, leg spread
- aero area, Cd and centre of pressure
- ascent and descent limits

**Two more limits:**
- The sim uses only the axial CG. The CAD's lateral offset (−42 mm, 6 mm) is kept in `IN.ctv.cg_lateral` but not flown.
- Monte Carlo scenarios (`mc_sim_setup`) model only the hopper so far, so `hopperWorkspace` rejects a CTV scenario.

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

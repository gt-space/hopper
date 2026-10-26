# Variants: one model, many configurations

`models/hopper_6dof.slx` is the variant harness. It replaces the full model copies that each added one feature (`_ground`, `_fluid_integration`, ...). Every place where those copies differ is a **Variant Subsystem**, and one workspace variable per axis picks the active choice. With every axis at its default, the harness is `hopper_6dof_NED_v2`, bit for bit.

## Running a configuration

```matlab
V = variantConfig('PROPULSION', 'fluid_network', 'GROUND', 'strut');
out = sim(buildSimInput('hopper_6dof', [], [], V));
```

`variantConfig` with no arguments gives the defaults. A misspelled axis or choice is an error, so a typo never silently runs the default. `buildSimInput` puts the `VARIANT_*` variables on the `SimulationInput`, and finds the variables the active choices read, so nothing goes in the base workspace.

## Axes

| Axis | Block | Choices (first is the default) |
|---|---|---|
| `SLOSH` | `Slosh` | `full`: slosh pendulums, `subsystems/slosh/Slosh.slx`. `off`: zero slosh force, moment and CG offsets. |
| `PROPULSION` | `Propulsion` | `lut`: thrust command passes through; `Thrust_to_Mdot_LUT` gives ox and fuel mdot. `fluid_network`: #26's Simscape feed system and engine give thrust and mdot. |
| `GROUND` | `Ground_Contact` | `none`: no ground reaction, as in v2. `strut`: #27's strut spring/damper per leg (`subsystems/ground/Ground_Contact.slx`), into the body force and moment sums. |
| `MASS_PROPS` | `Propellant_Flow` | `varying`: propellant drains, so mass, CG, MOI and dI/dt change. `constant`: zero propellant flow, so mass, CG and MOI hold their initial values and dI/dt = 0. Thrust is unaffected. |

Combinations to know:
- **`fluid_network` needs `GROUND=strut`.** The feed system takes about a second to build thrust, and without a ground the vehicle falls through the pad first.
- **`MASS_PROPS=constant`** is for vehicles that burn no propellant (the CTV). Propellant depletion cutoff then never triggers. The CG still moves sideways by under 1 mm as the engine gimbals, because the engine's own mass moves with the TVC.

Still to come: `VEHICLE` (hopper/CTV), `FDIR`, `FAULT_INJECTION`, `NAVIGATION`, `CONTROLLER`, `ENV`, `TANK_GEOM`.

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

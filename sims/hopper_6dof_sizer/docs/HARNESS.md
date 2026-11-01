# The harness: `models/hopper_6dof.slx`

The one top-level model. Every vehicle, physics option and GNC configuration is a variant choice inside it, not a copy of it. This page describes the model's structure. [VARIANTS.md](VARIANTS.md) covers the choices, and the [User guide](USER_GUIDE.md) covers running it.

## Top level: seven groups, six buses

Each group is a plain (virtual) subsystem holding the original v2 blocks unchanged. Grouping them did not change the simulation: the default flight is bit-identical to before. Each group sends one bus, and each bus element carries a short snake_case name shown at both ends (`dyn_bus . theta` where it leaves Vehicle_Dynamics, and again where Logging picks it up).

| Group | What's inside | Sends | Reads |
|---|---|---|---|
| `Environment` | `Wind_Forces` (shared file), the body-axis rotations of wind force and moment, the random roll-torque disturbance (`Random Source` × `Bernoulli Binary Generator`) | `env_bus` | `dyn_bus`, `prop_bus` |
| `Navigation_Sensors` | `Sensors` (shared sensor model; active when FDIR, fault injection or navigation is on), `Navigation` variant (`None` passes truth; `GMEKF_Shadow` / `GMEKF` run the nav team's filter), the truth-state mux | `nav_bus` | `dyn_bus`, `prop_bus`, `act_bus`, `gnc_bus` |
| `GNC` | `Controller` (CONTROLLER variant: LQR scheduler, LQI, legacy), `FDIR` variant, `Actuator_Faults` (fault injection on commands) | `gnc_bus` | `nav_bus`, `dyn_bus`, `prop_bus`, `act_bus` |
| `Actuation` | command demux, the thrust command chain (IC, rate limit, saturation, delay, quantizer), `Actuators_TVC_RCS` (shared file), weight from gravity | `act_bus` | `gnc_bus`, `prop_bus`, `dyn_bus` |
| `Mass_Propulsion` | `Landing_Engine_Cutoff`, `Propellant_Depletion_Cutoff`, `Propulsion` (PROPULSION variant), `Propellant_Flow`, `Propellant_Mass_Integration`, `Mass_Properties` (MASS_PROPS variant), total mass | `prop_bus` | `act_bus`, `dyn_bus` |
| `Vehicle_Dynamics` | `6DOF (Quaternion)`, `Slosh` (SLOSH variant), `Ground_Contact` (GROUND variant), force and moment sums, `Touchdown_Stop`, the Euler-angle and position taps | `dyn_bus` | `prop_bus`, `act_bus`, `env_bus` |
| `Logging` | To Workspace logs, nine `Scopes_<group>` variants, `Signal_Dump` | — | all six |

## Bus elements

Every signal that crosses between groups. The names are what you type in `scopeGroups.m` and what the signal log holds (`d.logsout.get('dyn_bus').Values.theta`).

**`env_bus`** (Environment)

| Element | Meaning |
|---|---|
| `wind_force` | wind force, body axes (N) |
| `wind_moment` | wind moment, body axes (N·m) |
| `roll_dist` | random roll-torque disturbance (N·m) |

**`gnc_bus`** (GNC)

| Element | Meaning |
|---|---|
| `u_cmd` | actuator command vector [thrust, TVC pitch, TVC yaw, RCS] after fault injection |
| `nav_mask` | FDIR's sensor mask to the nav filter |
| `nav_r_scale` | FDIR's measurement-noise inflation to the nav filter |

**`act_bus`** (Actuation)

| Element | Meaning |
|---|---|
| `force` | actuator force, body axes (N) |
| `moment` | actuator moment, body axes (N·m) |
| `tvc_pitch`, `tvc_yaw` | achieved gimbal angles |
| `thrust_cmd` | thrust command after rate limit, saturation, delay and quantization (N) |
| `thrust_req` | thrust requested by the controller (N) |

**`prop_bus`** (Mass_Propulsion)

| Element | Meaning |
|---|---|
| `thrust` | engine thrust (N) |
| `ox_mass`, `fu_mass` | remaining oxidizer and fuel (kg) |
| `mass` | vehicle mass (kg) |
| `cg` | CG, body axes from the base reference point (m) |
| `moi`, `moi_dot` | inertia tensor and its rate (3×3) |
| `engine_cg` | engine CG, which moves with the gimbal (m) |
| `cg_slosh`, `moi_slosh` | CG and inertia including the slosh masses |

**`nav_bus`** (Navigation_Sensors)

| Element | Meaning |
|---|---|
| `meas`, `valid`, `stamp` | the shared sensor model's measurement vector, validity flags and timestamps |
| `x_nav` | the 13-state the controller flies on: [pos NED; vel NED; body rates; quaternion]. Truth, or the estimate when `NAVIGATION=gmekf`. |
| `dcm_nav` | NED → body DCM matching `x_nav` |
| `x_filt`, `p_cov` | the GMEKF's 20-state and 19×19 covariance |
| `nav_diag` | per-sensor accept flags and NIS out of the correct step |
| `shadow` | the shadow vertical estimator |
| `x_vert` | a 13-state with the shadow estimator's vertical channel, which the R3 nav-degraded abort flies on |

**`dyn_bus`** (Vehicle_Dynamics)

| Element | Meaning |
|---|---|
| `vel_ned`, `pos_ned` | truth velocity and position, NED |
| `quat`, `dcm_be`, `w_body` | attitude quaternion, NED → body DCM, body rates |
| `x`, `y`, `z` | position components as logged (`out.x`, `out.y`, `out.z`; `z` is altitude) |
| `phi`, `theta`, `psi` | Euler angles (rad) |
| `u`, `w`, `p` | body velocity components and roll rate |
| `f_body` | specific force, body axes, including slosh |
| `cg_accel` | CG acceleration |
| `g_body` | gravity, body axes |
| `ox_slosh_pos`, `ox_slosh_mass`, `fu_slosh_pos`, `fu_slosh_mass` | slosh slug positions and masses |

## The Logging block

The left column holds the **To Workspace logs**: `out.x`, `out.y`, `out.z`, `out.velocity`, `out.u`, `out.roll`, `out.pitch`, `out.yaw`, `out.thrust`, `out.cg`, `out.cg_with_slosh`, `out.MOI_with_slosh`, `out.ox_mass`, `out.fuel_mass`.
- They are always written.
- Tests, Monte Carlo (`mc_main`) and analysis scripts read them, so don't rename them.

The right column holds the **switchable parts**. Each whole bus enters through an In Bus Element with no element selected (`env_bus_all`, ...) and feeds:
- **`Scopes_<group>`:** nine Variant Subsystems, each with one choice `On`, active when `VARIANT_SCOPES_<GROUP> == 2`. Inside: Bus Selectors, any conversions (deg, matrix diagonal, element slices), and one multi-display scope that opens when the simulation starts. The groups are defined in `src/setup/scopeGroups.m`.
- **`Signal_Dump`:** one choice `On`, active when `VARIANT_LOG == 2`. Logs each of the six buses whole to `logsout`.

`variantConfig('SCOPES', ..., 'LOG', ...)` sets these, and `buildSimInput` streams `logsout` to a MAT file when the dump is on. A choice that's off isn't compiled.

## Variants inside the groups

Variant activation time is *update diagram*: only the active choice is compiled, and only its variables are needed.

| Axis | Block |
|---|---|
| VEHICLE | none: a parameter set chosen in setup (`hopperWorkspace`) |
| SLOSH | `Vehicle_Dynamics/Slosh` |
| GROUND | `Vehicle_Dynamics/Ground_Contact` |
| PROPULSION | `Mass_Propulsion/Propulsion` |
| MASS_PROPS | `Mass_Propulsion/Mass_Properties` (and the propellant flow) |
| CONTROLLER | `GNC/Controller` |
| FDIR | `GNC/FDIR` |
| FAULT_INJECTION | `GNC/Actuator_Faults` and the sensor gate in `Navigation_Sensors/Sensors` |
| NAVIGATION | `Navigation_Sensors/Navigation` |
| SCOPES, LOG | `Logging/Scopes_*`, `Logging/Signal_Dump` |

## Changing the harness

**Small edits** (a gain, a block parameter, a new block inside a group): edit the model in the team release and save. Prefer editing the shared subsystem file under `subsystems/` when the change belongs there.

**Structural or repeated edits:** write a rerunnable script in `tools/refactor/` and run it on the model, so the change can be reviewed and redone on a fresh copy. The scripts so far:

| Script | What it did |
|---|---|
| `regroupHarness` | Deleted 44 commented-out leftovers. Moved every top-level block into the seven groups (`moveIntoGroup`). Converted the 76 wires between groups into 21 bus lines with named elements (`busifyHarness`). Laid out the result (`layoutHarnessGroups`). |
| `moveIntoGroup` | Moves blocks into a new subsystem. It records every connection, then redraws it through named ports, never deleting a live line. (`Simulink.BlockDiagram.createSubsystem` changed the flight and crashed MATLAB, so it isn't used.) |
| `busifyHarness` | One bus per group. The element names come from the table in its `signalNames()`. It asserts that every group output has a name and that the names are unique. |
| `layoutHarnessGroups` | `arrangeSystem` inside each group and, once on buses, at the top level. |
| `buildLoggingScopes` | Builds the scope groups and the signal dump from `scopeGroups.m`. Rerun it after editing that list. |
| `fixTvcSkewAxes` | The engine-CG swing in `Mass_Properties` used yaw and pitch swapped relative to the thrust model. Fixed, with a wiring check. |
| `parameterizeVehicleGeometry`, `parameterizeWindForces`, `parameterizeLandingCutoff`, `parameterizeRollActuation` | Replaced literals in shared subsystems with vehicle-config and mission values. |
| `fixBernoulliSeeds` | Gave each gust generator its own fixed seed, for repeatable runs. |
| `updateNavigationFilters`, `syncNavigation` | The nav team's filter changes the harness needs, and merging their newer work. |

**Adding a signal between groups:**
1. Wire it to an Out Bus Element on the sending group's bus. Copy an existing one inside that group and change its `Element`; a new block can't take an existing port name.
2. Pick it up with an In Bus Element where it's needed.
3. Name it in snake_case, short but unambiguous, and add it to the tables above.

**Checking a structural change:** compare every logged output of the default flight before and after, as `tools/refactor` changes so far have done (bit-identical), then run `buildtool check test sim`.

**Block paths changed with the regroup.** Code that reaches into the harness by path must use the group path:
- `[model '/Mass_Propulsion/Mass_Properties']`, not `[model '/Mass_Properties']`
- `massPropertiesCheck`, `sloshOpenLoopTest` and `parameterizeRollActuation` are already updated.

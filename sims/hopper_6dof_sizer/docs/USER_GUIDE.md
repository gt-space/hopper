# User guide

How to run, configure, debug and analyse the hopper 6DOF sim. The [README](../README.md) has the short version. [HARNESS.md](HARNESS.md) describes the model, and [KNOWN_ISSUES.md](KNOWN_ISSUES.md) lists what's broken or unfinished.

- [Setup](#setup)
- [Running a flight](#running-a-flight)
- [Reading the results](#reading-the-results)
- [Configurations (variants)](#configurations-variants)
- [Changing a parameter for one run](#changing-a-parameter-for-one-run)
- [Scopes and the signal log](#scopes-and-the-signal-log)
- [Running from the Simulink window](#running-from-the-simulink-window)
- [Monte Carlo](#monte-carlo)
- [FDIR and fault injection](#fdir-and-fault-injection)
- [Navigation](#navigation)
- [The CTV](#the-ctv)
- [Changing the vehicle or the mission](#changing-the-vehicle-or-the-mission)
- [Recipes](#recipes)
- [Troubleshooting](#troubleshooting)

## Setup

Do this once per machine:

1. **Install MATLAB R2026a** (the release in [`TEAM_RELEASE`](../TEAM_RELEASE)) with these products:
   - Simulink, Aerospace Blockset
   - DSP System, Communications, Signal Processing, Control System, Optimization, and Statistics and Machine Learning Toolbox
   - Simscape and Simscape Fluids, for `PROPULSION=fluid_network`
   - Parallel Computing Toolbox, optional, for parallel Monte Carlo

   Other releases can sit alongside, but only **save** models from R2026a; CI rejects anything newer. If `matlab` on your PATH is a newer release, start R2026a from its own shortcut or its full path, for example `C:\Program Files\MATLAB\R2026a\bin\matlab.exe`.
2. **Install Python with CoolProp**, needed only for Monte Carlo (`prop_system` calls it every run). In MATLAB, `pyenv` should show your Python.
   ```bash
   pip install CoolProp
   ```
3. **Set up the git merge tools for `.slx` files:** the three `git config` lines in the README.

Every session: open `HopperSim.prj`, by double-clicking it or with `openProject('sims/hopper_6dof_sizer')`. The project:
- puts `src/`, `subsystems/`, `inputs/` and the other folders on the path
- sends `slprj` and codegen caches to `work/`
- starts Python if one is configured

Don't `cd` into folders or `addpath`.

## Running a flight

```matlab
out = sim(buildSimInput());                                   % nominal hopper
out = sim(buildSimInput('hopper_6dof', [], [], V));           % a configuration, V from variantConfig
out = sim(buildSimInput('hopper_6dof', scenario));            % one Monte Carlo scenario
out = sim(buildSimInput('hopper_6dof', [], W, V));            % your own workspace struct W
```

**What `buildSimInput(model, scenario, W, variants)` does:**
1. Builds the workspace struct `W` with `hopperWorkspace(scenario, variants)`:
   - the vehicle config and mission inputs
   - mass-property tables
   - the reference trajectory and LQR gains (`LinerizationMaster`)
   - the baffle-damping and mdot lookup tables
2. Adds the sensor, navigation and FDIR parameters when those variants are on.
3. Finds which variables the model reads (`modelVariables`, cached per variant combination) and sets each one on a `Simulink.SimulationInput`.

The run needs nothing in your base workspace, and two inputs can run side by side with `parsim`.

**Timing:**
- The first `sim` in a session compiles the model (a minute or two).
- Building `W` takes a few seconds. To fly several runs on the same setup, build it once with `W = hopperWorkspace([], V)` and pass it in.

**Stop time:** the flight ends at touchdown (`Touchdown_Stop`) or at the end of the trajectory. With FDIR on, the stop time is the trajectory end + 30 s, so an abort can finish. For a short run:
```matlab
out = sim(buildSimInput().setModelParameter('StopTime', '5'));
```

## Reading the results

`out` is a `Simulink.SimulationOutput`. The harness always writes these timeseries (`out.<name>.Time`, `.Data`):

| Variable | Signal |
|---|---|
| `x`, `y`, `z` | position (m); `z` is altitude |
| `velocity` | NED velocity (m/s) |
| `u` | body velocity component |
| `roll`, `pitch`, `yaw` | Euler angles: `roll` and `pitch` in deg, `yaw` in rad (as v2 logged them) |
| `thrust` | engine thrust (N) |
| `ox_mass`, `fuel_mass` | remaining propellant (kg) |
| `cg`, `cg_with_slosh`, `MOI_with_slosh` | mass properties |
| `tout` | time vector |

Some variants add their own logs:
- **FDIR:** `out.get('fdir_mode')`, `out.get('fdir_level')`, ...
- **Navigation:** `out.get('nav_estimate')`, `out.get('nav_error')`, `out.get('nav_cov_diag')`

`out.who` lists everything that was logged.

```matlab
plot(out.x.Data, out.y.Data), axis equal, title('ground track')
fprintf('landing offset %.3f m\n', hypot(out.x.Data(end), out.y.Data(end)))
```

For anything not in that list, use the [signal log](#scopes-and-the-signal-log).

## Configurations (variants)

`variantConfig` turns names into the control variables the model checks. Pass it to `buildSimInput`:

```matlab
V = variantConfig('SLOSH', 'off', 'GROUND', 'strut');
V = variantConfig(struct('VEHICLE', 'ctv', 'NAVIGATION', 'gmekf'));
```

| Axis | Choices (first = default) |
|---|---|
| `VEHICLE` | `hopper`, `ctv` |
| `SLOSH` | `full`, `off` |
| `PROPULSION` | `lut`, `fluid_network` (Simscape feed system) |
| `GROUND` | `none`, `strut` |
| `MASS_PROPS` | `varying`, `constant` |
| `CONTROLLER` | `lqr_sched`, `lqi` (experimental), `lqr_sched_legacy` |
| `FDIR` | `off`, `on` |
| `FAULT_INJECTION` | `off`, `on` |
| `NAVIGATION` | `none` (truth), `gmekf_shadow`, `gmekf` |

Logging options sit beside the axes: `SCOPES` and `LOG` (see below).

- **Errors on typos:** an unknown axis or choice is an error, listing the valid ones.
- **CTV defaults:** `VEHICLE='ctv'` sets `SLOSH='off'` and `MASS_PROPS='constant'` unless you name them.
- **More detail:** [VARIANTS.md](VARIANTS.md) explains each choice.

## Changing a parameter for one run

Build `W` yourself, change it, and pass it in:

```matlab
V = variantConfig();
W = hopperWorkspace([], V);
W.IN.wind.uwind = 8;  W.IN.wind.vwind = 0;     % steady 8 m/s north wind
out = sim(buildSimInput('hopper_6dof', [], W, V));
```

**What this reaches:** it works for anything the model reads directly. **It doesn't** recompute what the setup derived from that value, such as the reference trajectory, LQR gains and mass-property tables. To change those, change the source and rebuild `W`:
- the vehicle config (`vehicles/<vehicle>/<vehicle>_config.m`)
- `inputs/mission_inputs.m`
- a Monte Carlo scenario

`runMonteCarlo(..., Variables=struct(...))` does the same per run.

## Scopes and the signal log

**Scope groups:** open live scopes while the sim runs.
```matlab
V = variantConfig('SCOPES', ["trajectory" "attitude" "actuation"]);
V = variantConfig('SCOPES', 'all');
```

| Group | Displays |
|---|---|
| `trajectory` | x and y; altitude vs `z_ref`; vertical speed vs `vz_ref`; NED velocity |
| `attitude` | roll, pitch and yaw (deg); body rates (deg/s); quaternion |
| `propulsion` | thrust: actual, command, request and profile; ox and fuel mass; vehicle mass |
| `actuation` | TVC pitch and yaw; body force; body moment; controller command `u` |
| `mass_props` | CG; CG with slosh; engine CG; MOI and dI/dt diagonals |
| `slosh` | slug positions; slug masses; body specific force; CG acceleration |
| `environment` | wind force; wind moment; roll disturbance |
| `navigation` | position and velocity, estimate vs truth; filter diagnostics; shadow estimator |
| `fdir` | nav sensor mask; nav R scale; sensor valid flags |

- **When they open:** each group's scope opens when the simulation starts. In `matlab -batch` or on parallel workers there is no display, and nothing opens.
- **Cost when off:** none. A group that's off isn't compiled, and all groups default to off.
- **Combining:** they combine with any axis, for example `variantConfig('NAVIGATION', 'gmekf', 'SCOPES', 'navigation')`.
- **Navigation estimate vs truth:** this is only meaningful with `NAVIGATION` set to `gmekf_shadow` or `gmekf`. With `none`, both lines are truth.

**Adding or changing a group:**
1. Edit `src/setup/scopeGroups.m`. Each display is a title and a list of signals:
   - `'<bus>.<element>'`, optionally followed by `:deg`, `:diag` or `:i1-3`
   - `'ref.<variable>'` for a reference timeseries
   - The element names are in [HARNESS.md](HARNESS.md#bus-elements).
2. Run `tools/refactor/buildLoggingScopes` in R2026a. It rebuilds every group.
3. Save the model and commit both files.

**Signal log:** every signal between the groups, written to a MAT file.
```matlab
V = variantConfig('LOG', 'on');                            % outputs/logs/hopper_<yyyyMMdd_HHmmss>.mat
V = variantConfig('LOG', 'outputs/logs/wind8.mat');         % or a file you name
out = sim(buildSimInput('hopper_6dof', [], [], V));

d = load(V.LOG_FILE);
theta = d.logsout.get('dyn_bus').Values.theta;              % a timeseries
plot(theta.Time, rad2deg(squeeze(theta.Data)))
Simulink.sdi.view                                          % or browse it in the Data Inspector
```

- **Contents:** the log holds the six buses as structs of timeseries, all 53 signals, at the solver step.
- **Size:** a nominal hopper flight is ~42 MB, mostly the nav covariance `nav_bus.p_cov`. Leave the log off in large Monte Carlos.
- **Workspace logs are unaffected:** `out.x` and the other workspace logs are written either way.

## Running from the Simulink window

The model reads its variables from the `SimulationInput`, not from the base workspace, so pressing Run on a freshly opened model fails with missing variables. To use the Simulink window, apply an input to the model first:

```matlab
in = buildSimInput('hopper_6dof', [], [], variantConfig('SCOPES', 'attitude'));
open_system('hopper_6dof');
in.applyToModel();      % the model now carries those variables and variant choices
% press Run (or sim('hopper_6dof')) as often as you like
```

When you're done, **close the model without saving**, since the applied state shouldn't be committed:
```matlab
close_system('hopper_6dof', 0)
```
`buildSimInput` doesn't work on a model with an input applied. Close it first.

**Older route:** `sim_setup`, which fills the base workspace, still exists for the `fsw_bridge` scripts. It predates the variants, so it doesn't set the variant controls the harness needs.

**Just looking:** to open the model without running it, use `open_system('hopper_6dof')`.

## Monte Carlo

```matlab
R = runMonteCarlo(10);                                     % 10 seeded runs + the nominal
R = runMonteCarlo(200, UseParallel=true);                  % parsim; pool sized to free RAM
R = runMonteCarlo(50, Variants=variantConfig('VEHICLE', 'ctv'));
R = runMonteCarlo(100, Wind=false, OutputFile=hopperRoot('results', 'mc_nowind.mat'));
```

**Scenarios:** `generateScenarios` draws them with Latin hypercube sampling from the bounds in:
- `inputs/mc_params.json` for the hopper
- `inputs/mc_params_ctv.json` for the CTV

The rng is seeded with 42, so run *k* is the same scenario in every session and on every machine.

**What you get back:**
- `R.summary`: one row per run, with apex, landing velocity, lateral velocity, landing offset, tilt and flight time.
  - `pass` is the result against the CI regression limits (`tests/montecarlo/mcLimits.m`).
  - `meets_requirements` and `missed_requirements` are the result against the mission requirements, the `[req]` fields in `inputs/mission_inputs.m`.
- `R.results`: full per-run results. `R.scenarios`: the scenarios. `R.nominal`: the nominal run.

**Options:** `Variants`, `Model`, `Wind`, `ParamsFile`, `Variables` (struct of model variables set after setup), `UseParallel`, `Workers`, `Nominal`, `Limits`, `OutputFile`. See `help runMonteCarlo`.

**Memory:** each parallel worker commits ~3.3 GB. `mcWorkerCount` picks a pool that fits, and `Workers` overrides it.

**Large campaigns:** thousands of runs, with checkpoints, use `mc_parallelization`. It is built from the same pieces.

### Debugging a Monte Carlo run

Every run can be reproduced exactly. To re-fly run *k* with every scope open and its signals logged:

```matlab
R = runMonteCarlo(10);
k = find(~R.summary.meets_requirements, 1);
V = variantConfig('SCOPES', 'all', 'LOG', sprintf('outputs/logs/mc_run%d.mat', k));
out = sim(buildSimInput('hopper_6dof', R.scenarios(k), [], V));
R.scenarios(k)          % the dispersions that made it fail
```

For the CTV, pass the same `VEHICLE` choice in `V` that the campaign used.

## FDIR and fault injection

```matlab
V = variantConfig('FDIR', 'on', 'FAULT_INJECTION', 'on');
W = hopperWorkspace([], V);
W.faults = faultScenario({'lidar', 'loss', 20, inf, 0, 2}, {'gimbal_pitch', 'stuck', 8});
out = sim(buildSimInput('hopper_6dof', [], W, V));
plot(out.get('fdir_mode'))      % PDR mode sequence
plot(out.get('fdir_level'))     % response level: -1 none, 0 R0 mask, 1 R1, ...
```

- **Fault rows:** each is `{target, mode, t_on, t_off, magnitude, param}`; `help faultScenario` has the targets and modes.
- **FDIR without injection:** `FDIR=on` alone flies the nominal mission through the mode manager. A quiet flight walks the PDR mode sequence to LANDED and never raises a response.
- **Injection without FDIR:** `FAULT_INJECTION=on` with `FDIR=off` shows the unmitigated response to a fault.
- **More detail:** [FDIR_ARCHITECTURE.md](FDIR_ARCHITECTURE.md) has the catalog (`inputs/fdir_fault_catalog.csv`), monitors, thresholds and responses. Add `'SCOPES', 'fdir'` to watch the masks.

## Navigation

```matlab
V = variantConfig('NAVIGATION', 'gmekf_shadow', 'SCOPES', 'navigation');   % filter runs, truth flies
V = variantConfig('NAVIGATION', 'gmekf');                                   % fly on the estimate
out = sim(buildSimInput('hopper_6dof', [], [], V));
plot(out.get('nav_error'))      % estimate - truth, in the controller's 13-state
```

- **Whose filter:** it's the nav team's (`subsystems/navigation`); our adapter is in `src/navigation`.
- **What feeds it:** the shared `Sensors` model, so FDIR and the filter see the same measurements.
- **Tuning:** to try other values, set `W.nav_Qc` or `W.nav_P0` before `buildSimInput`.
- **More detail:** [VARIANTS: Navigation](VARIANTS.md#navigation).

## The CTV

```matlab
out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('VEHICLE', 'ctv')));
R = runMonteCarlo(10, Variants=variantConfig('VEHICLE', 'ctv'));
```

- **Vehicle and mission:** the vehicle is `vehicles/ctv/ctv_config.m`; its mission and controller weights are `inputs/ctv_inputs.m`.
- **Nominal flight:** in calm air it meets every requirement and lands on the estimate too.
- **Dispersed runs:** mostly fail (see [KNOWN_ISSUES.md](KNOWN_ISSUES.md)).
- **`ASSUMPTION` values:** many CTV values are marked `ASSUMPTION` and need confirming before a CTV result means much.
- **More detail:** [VARIANTS: The CTV](VARIANTS.md#the-ctv).

## Changing the vehicle or the mission

| Change | Where |
|---|---|
| Masses, CGs, inertia, tanks, engine, legs, RCS, sensors, aero | `vehicles/<vehicle>/<vehicle>_config.m` ([VEHICLE_CONFIG.md](VEHICLE_CONFIG.md)) |
| Requirements (`[req]` fields), trajectory, wind, cutoff height | `inputs/mission_inputs.m` (CTV overrides in `inputs/ctv_inputs.m`) |
| Monte Carlo dispersions | `inputs/mc_params.json`, `inputs/mc_params_ctv.json` |
| CI regression limits | `tests/montecarlo/mcLimits.m` (owned by the GNC leads) |
| FDIR catalog | `inputs/fdir_fault_catalog.csv` |
| LQR weights | `IN.gnc.lqr` in the mission or CTV inputs |

**Slosh baffles:** after changing tanks or baffles, rerun `generateBaffleLUT`. The setup warns when the LUT is stale.

**Flight changes on purpose:** if a change is meant to move the nominal flight, run `NominalFlightTest.updateBaseline()` and say so in the PR.

## Recipes

**Compare two configurations:**
```matlab
W = hopperWorkspace();
a = sim(buildSimInput('hopper_6dof', [], W, variantConfig()));
b = sim(buildSimInput('hopper_6dof', [], W, variantConfig('SLOSH', 'off')));
plot(a.z.Time, a.z.Data, b.z.Time, b.z.Data), legend('slosh', 'no slosh')
```
Sharing `W` saves the setup time. Don't share it across `VEHICLE` choices; each vehicle needs its own.

**Sweep a parameter:**
```matlab
V = variantConfig(); W0 = hopperWorkspace([], V);
winds = 0:2:10;
ins = Simulink.SimulationInput.empty;
for k = 1:numel(winds)
    W = W0; W.IN.wind.uwind = winds(k);
    ins(k) = buildSimInput('hopper_6dof', [], W, V);
end
outs = parsim(ins, 'TransferBaseWorkspaceVariables', 'off');   % or sim(ins) serially
offset = arrayfun(@(o) hypot(o.x.Data(end), o.y.Data(end)), outs);
```

**Check before a pull request:** `buildtool check test` (~2 min) on every change. Add `buildtool sim` (~40 min) when the change touches the model, a subsystem, setup or GNC code, and `buildtool mc` when it might move the dispersed results. The PR template asks which you ran.

**Find where a signal lives:** [HARNESS.md](HARNESS.md#bus-elements) lists every bus element. Inside the model, the In Bus Element blocks show `bus . element`.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `'X' reads 'Y', which the setup did not create` | A model variable that no setup step makes, often from a new block parameter. Create it in the vehicle config, the mission inputs or `hopperWorkspace`. |
| Run button: undefined variable `VARIANT_...` or `IN` | The model isn't fed from the base workspace. See [Running from the Simulink window](#running-from-the-simulink-window). |
| `buildSimInput` fails in `Simulink.findVars` | The model has an input applied (`applyToModel`) or unsaved edits that break compilation. `close_system('hopper_6dof', 0)` and retry. |
| `y = logical ...` printed many times | Known console noise, harmless ([KNOWN_ISSUES.md](KNOWN_ISSUES.md)). |
| CI `guards` fails: model newer than the team release | Saved from a newer MATLAB. Re-save with `Simulink.exportToVersion(model, 'tmp.slx', teamRelease())` and replace the file. |
| A `.slx.r2026a` file appeared | Simulink backed up a model it upgraded on save, meaning you saved from a newer release. Delete the backup and re-save from R2026a. |
| `Unable to create a test suite` | Run `openProject` first: the tests need the project path. |
| Monte Carlo fails with a Python or CoolProp error | Install CoolProp into the Python `pyenv` points at. |
| Parallel pool runs out of memory | Pass `Workers=` with a smaller number. Each worker needs ~3.3 GB. |
| A file you added isn't found on another machine | Run `syncProjectFiles` so the project lists it, and commit `resources/project`. |
| "Model was exported from R2026b" message in tests | `throttle_valve_simulink.slx` was saved in R2026b ([KNOWN_ISSUES.md](KNOWN_ISSUES.md)). Harmless; don't save it from R2026a by accident. |

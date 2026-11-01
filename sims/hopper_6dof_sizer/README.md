# Hopper 6DOF Simulation

MATLAB/Simulink 6DOF simulation of the Hopper (YJSP's propulsive lander) and its Controls Test Vehicle (CTV).

It models:
- **Vehicle:** rigid-body dynamics, propellant depletion with varying mass, CG and inertia, and propellant slosh.
- **Actuation:** TVC and RCS, the throttle valve, and optionally a Simscape feed system.
- **Ground and environment:** ground contact and wind.
- **GNC:** the LQR controller, the nav team's GMEKF navigation filter, and FDIR with fault injection.

It also runs seeded Monte Carlo campaigns against the mission requirements.

> **New here?** Read the [Quick start](#quick-start), then the [User guide](docs/USER_GUIDE.md).
> **Changing the sim?** Read [CONTRIBUTING](docs/CONTRIBUTING.md) first.
> **Something looks wrong?** Check [Known issues](docs/KNOWN_ISSUES.md).

## Quick start

1. **Use MATLAB R2026a.** The release is named in [`TEAM_RELEASE`](TEAM_RELEASE). CI runs it, and models must be saved in it.
   - If you have several releases, check which one opens. On a machine where `matlab` on the PATH is a newer release, start R2026a by its full path.
2. **Open the project.** Double-click `HopperSim.prj`, or run:
   ```matlab
   openProject('sims/hopper_6dof_sizer')
   ```
   This sets the whole path and sends Simulink caches to `work/`. No `cd` or `addpath` is needed.
3. **Fly the nominal hopper:**
   ```matlab
   out = sim(buildSimInput());
   plot(out.z.Time, out.z.Data), xlabel('t (s)'), ylabel('altitude (m)')
   ```
4. **Fly another configuration, with scopes open:**
   ```matlab
   V = variantConfig('VEHICLE', 'ctv', 'SCOPES', ["trajectory" "attitude"]);
   out = sim(buildSimInput('hopper_6dof', [], [], V));
   ```
5. **Run a Monte Carlo.** This needs Python with CoolProp: `pip install CoolProp`.
   ```matlab
   R = runMonteCarlo(10);
   R.summary
   ```
6. **Before you push,** run what CI runs:
   ```matlab
   buildtool check test sim
   ```

The first simulation in a session compiles the model and takes a minute or two; later ones reuse the cache. The [User guide](docs/USER_GUIDE.md) covers every option.

## What the refactor changed

The sim used to be a flat folder of ~120 files. Every feature was a full copy of the top-level model (`_ground`, `_fault_injection`, `_fluid_integration`, ...), and the setup scripts were joined through the base workspace. There were no tests and no CI. Now:

| Area | Now |
|---|---|
| **One model** | `models/hopper_6dof.slx` is the only top-level model. Features are **variant choices** in it, not copies. |
| **Readable harness** | The top level is seven groups joined by six named buses. Every signal between groups has a short snake_case name, such as `dyn_bus.theta` or `prop_bus.cg`. See [HARNESS.md](docs/HARNESS.md). |
| **Shared subsystems** | Slosh, mass properties, actuators, wind, ground contact, propulsion, FDIR and navigation live in their own files under `subsystems/`. Edit one file and every user of it changes. |
| **One setup call** | `buildSimInput(model, scenario, W, variants)` returns a `Simulink.SimulationInput` carrying every variable the model reads. Nothing has to be in the base workspace, and `parsim` can run inputs side by side. |
| **Variant axes** | `variantConfig` selects VEHICLE, SLOSH, PROPULSION, GROUND, MASS_PROPS, CONTROLLER, FDIR, FAULT_INJECTION and NAVIGATION. A typo is an error, not a silent default. See [VARIANTS.md](docs/VARIANTS.md). |
| **Vehicle config files** | `vehicles/<vehicle>/<vehicle>_config.m` is the one source of geometry, masses and hardware. The conceptual sizer no longer runs in setup. See [VEHICLE_CONFIG.md](docs/VEHICLE_CONFIG.md). |
| **CTV as a vehicle** | `variantConfig('VEHICLE', 'ctv')`: its own CAD mass properties, EDF thrust, roll props and mission, flown in the same harness. |
| **Navigation** | The nav team's GMEKF runs in the harness, fed by the shared sensor model. `gmekf_shadow` runs it beside truth; `gmekf` flies on it. The hopper and the CTV both land on the estimate. |
| **FDIR and fault injection** | The PDR mode manager, sensor, GNC and estimator monitors, responses R0–R4, and a fault table that drives sensor and actuator faults. See [FDIR_ARCHITECTURE.md](docs/FDIR_ARCHITECTURE.md). |
| **Debug scopes and logs** | Nine scope groups you can switch on, and a log of every bus signal to a MAT file: `variantConfig('SCOPES', ..., 'LOG', 'on')`. Both are off by default and cost nothing then. |
| **Monte Carlo** | `runMonteCarlo(n, ...)`: seeded, serial or `parsim`, scored against CI regression limits and, separately, against the mission requirements. |
| **Tests and CI** | `buildtool check test sim mc` runs locally and on every push and pull request. CI also produces a model-diff report for changed `.slx` files. |
| **Reproducible harness edits** | Every structural change to the harness is a rerunnable script in `tools/refactor/`, for example `regroupHarness`, `buildLoggingScopes` and `fixTvcSkewAxes`. |
| **Team process** | CODEOWNERS teams, a pull request template asking for variant and Monte Carlo results, a team MATLAB release, and Simulink-aware diff and merge for `.slx`. |

## The harness at a glance

```
 Environment ──env_bus──┐
 Navigation_Sensors ─nav_bus─┐                    (every group reads the buses it needs;
 GNC ──gnc_bus──┐            │                     Logging reads all six)
 Actuation ──act_bus──┐      │
 Mass_Propulsion ──prop_bus──┤
 Vehicle_Dynamics ──dyn_bus──┘──► Logging (workspace logs, scope groups, signal dump)
```

| Group | Contains | Sends |
|---|---|---|
| `Environment` | wind forces and moments (`Wind_Forces`), random roll disturbance | `env_bus` |
| `Navigation_Sensors` | shared `Sensors` model, `Navigation` variant (truth or GMEKF) | `nav_bus` |
| `GNC` | `Controller` (LQR variants), `FDIR` variant, `Actuator_Faults` | `gnc_bus` |
| `Actuation` | thrust command chain, `Actuators_TVC_RCS` | `act_bus` |
| `Mass_Propulsion` | engine cutoffs, `Propulsion` variant, propellant integration, `Mass_Properties` | `prop_bus` |
| `Vehicle_Dynamics` | 6DOF, `Slosh` variant, `Ground_Contact` variant, force and moment sums | `dyn_bus` |
| `Logging` | To Workspace logs, `Scopes_<group>` variants, `Signal_Dump` | none |

[HARNESS.md](docs/HARNESS.md) has every bus element, where each group's blocks came from, and how to change the harness.

## Common tasks

| I want to... | Do this | More |
|---|---|---|
| Fly the nominal hopper | `out = sim(buildSimInput())` | [User guide](docs/USER_GUIDE.md#running-a-flight) |
| Fly the CTV | `variantConfig('VEHICLE', 'ctv')` | [VARIANTS: The CTV](docs/VARIANTS.md#the-ctv) |
| Turn off slosh, add ground struts, ... | `variantConfig('SLOSH', 'off', 'GROUND', 'strut')` | [VARIANTS](docs/VARIANTS.md#axes) |
| Watch signals while it flies | `variantConfig('SCOPES', ["attitude" "actuation"])` | [User guide](docs/USER_GUIDE.md#scopes-and-the-signal-log) |
| Save every signal for later | `variantConfig('LOG', 'on')` | [User guide](docs/USER_GUIDE.md#scopes-and-the-signal-log) |
| Fly on the navigation filter | `variantConfig('NAVIGATION', 'gmekf')` | [VARIANTS: Navigation](docs/VARIANTS.md#navigation) |
| Inject a fault and watch FDIR | `W.faults = faultScenario(...)`, with `FDIR` and `FAULT_INJECTION` on | [User guide](docs/USER_GUIDE.md#fdir-and-fault-injection) |
| Run a Monte Carlo | `R = runMonteCarlo(100, UseParallel=true)` | [User guide](docs/USER_GUIDE.md#monte-carlo) |
| Re-fly one Monte Carlo case with scopes | `sim(buildSimInput('hopper_6dof', R.scenarios(k), [], variantConfig('SCOPES', 'all')))` | [User guide](docs/USER_GUIDE.md#debugging-a-monte-carlo-run) |
| Change a mass, tank or leg | Edit `vehicles/hopper/hopper_config.m` | [VEHICLE_CONFIG](docs/VEHICLE_CONFIG.md) |
| Change a requirement | `inputs/mission_inputs.m`, fields marked `[req]` | [CONTRIBUTING](docs/CONTRIBUTING.md#requirements-and-ci-limits) |
| Add a feature or a variant choice | A subsystem file plus a variant choice | [VARIANTS](docs/VARIANTS.md#adding-a-choice-or-an-axis), [CONTRIBUTING](docs/CONTRIBUTING.md#adding-a-feature) |
| Press Run in the Simulink window | `buildSimInput(...).applyToModel()` first | [User guide](docs/USER_GUIDE.md#running-from-the-simulink-window) |

## Testing

From the project root (`buildtool` needs R2023b or later):

| Command | Runs | Time |
|---|---|---|
| `buildtool` | `check` + `test`: repo guards; unit, model-file and FDIR tests; the nav team's filter tests | ~2 min |
| `buildtool sim` | Nominal flight against `tests/baselines/nominal.mat`, the slosh open-loop check, a 2 s smoke run of every variant choice, the scopes and log check, and the FDIR and navigation flights | ~40 min |
| `buildtool mc` | Seeded 10-run hopper Monte Carlo against `tests/montecarlo/mcLimits.m`, plus 5 seeded CTV runs | ~5 min with a parallel pool |

- **What fails `check`:** committed build caches or autosaves, or a model saved in a release newer than `TEAM_RELEASE`.
- **What fails `mc`:** a run that errors, diverges or does worse than today's vehicle plus margin. Separately, `R.summary.meets_requirements` scores each run against the mission requirements, which are stricter; today most dispersed hopper runs miss the landing-offset requirement (see [Known issues](docs/KNOWN_ISSUES.md)).
- **When your change is meant to move the flight:** run `NominalFlightTest.updateBaseline()`, commit the new baseline and say so in the PR.
- **Repeatability:** the same inputs give bit-identical results in any session.

### CI

[`.github/workflows/hopper-sim.yml`](../../.github/workflows/hopper-sim.yml) runs on every push and pull request that touches this folder. It uses GitHub-hosted Linux runners with the release in `TEAM_RELEASE`; the repo is public, so MathWorks provides the license for these runners.

| Job | What it does |
|---|---|
| `guards` | No MATLAB, finishes in seconds. `tools/ci/guards.sh`: no committed caches, no model newer than the team release. `tools/ci/codeowners_check.py`: every owner in `.github/CODEOWNERS` resolves. |
| `test` | `buildtool check test sim`. Uploads JUnit XML as the `test-results` artifact. |
| `montecarlo` | `buildtool mc`. Writes the hopper and CTV per-run tables to the run's summary page and uploads `outputs/montecarlo/`. |
| `model-diff` | Pull requests only. A Simulink Comparison HTML report for every changed `.slx`, uploaded as `model-diffs`. |

To make `main` require these jobs, an admin turns this on once under **Settings → Branches → Add branch protection rule** for `main`:
- **Require status checks to pass:** `Repo guards`, `Unit, model and nominal-flight tests`, `Monte Carlo (10 seeded runs)`.
- **Require a pull request before merging,** with 1 approval.

## Layout

| Folder | Contents |
|---|---|
| `models/` | `hopper_6dof.slx`, the master variant harness. `hopper_6dof_NED_v2_fswBridge.slx`, the flight-software bridge. `dev/`: work in progress (LQI). `standalone/`: component models. |
| `subsystems/` | Shared subsystem files the harness references: `slosh` (`Slosh`, plus the `Slosh_Axis` pendulum used six times), `mass_properties`, `actuators` (TVC/RCS, EDF fans), `environment` (wind), `ground`, `gnc` (LQR gain scheduler), `propulsion`, `flight_logic` (cutoffs, touchdown stop), `fdir`, `navigation` (the nav team's filters). |
| `src/` | MATLAB code by area:<br>• setup and design: `setup` (`buildSimInput`, `variantConfig`, `scopeGroups`, ...), `sizing` (design trades; not run in setup)<br>• physics: `dynamics`, `propulsion`, `slosh`<br>• GNC and FDIR: `gnc`, `sensors` (the shared sensor model), `navigation` (the adapter to the nav filter), `fdir`, `faults`<br>• analysis: `monte_carlo`, `analysis`, `util` |
| `vehicles/` | `hopper/hopper_config.m`, `ctv/ctv_config.m`: geometry, masses, hardware. |
| `inputs/` | `mission_inputs.m` (requirements, trajectory, environment), `ctv_inputs.m`, `mc_params.json`, `mc_params_ctv.json`, `fdir_fault_catalog.csv`. |
| `config/` | `hopper_harness.sldd`, the harness data dictionary (references the nav team's bus types). |
| `data/` | Lookup tables: `mdot_lookup.xlsx`, `wind_vectors2.mat`, `cg_I_LUT.mat`, `baffle_lut.*`, `engine_contour.xlsx`. |
| `results/` | Tracked Monte Carlo and sweep results (`.mat`). |
| `tests/` | `unit`, `models`, `fdir`, `sim`, `montecarlo`, `baselines`. |
| `validation/` | Validation studies: mass properties, slosh open loop. |
| `fsw_bridge/` | C codegen and Rust SIL bridge to flight software. |
| `tools/refactor/` | Rerunnable scripts that change the harness or shared files. `tools/ci/`: CI helpers. |
| `docs/` | The documents listed below. |
| `archive/` | Superseded models and scripts, kept for reference and **not** on the path. |
| `outputs/`, `work/` | Generated figures, logs, CSVs and build caches. Gitignored. |

- **Moving files:** if you move or add files with `git mv` or Explorer instead of the Project panel, run `syncProjectFiles`.
- **Writing files:** scripts that write files use `hopperRoot('data'|'results'|'outputs', ...)`.

## Documents

| Document | Read it for |
|---|---|
| [USER_GUIDE.md](docs/USER_GUIDE.md) | How to run, configure, debug and analyse: every workflow, with examples |
| [HARNESS.md](docs/HARNESS.md) | The model's structure: groups, buses, every signal name, the Logging block, the refactor scripts |
| [VARIANTS.md](docs/VARIANTS.md) | Every variant axis and choice, navigation, the CTV, adding a choice |
| [VEHICLE_CONFIG.md](docs/VEHICLE_CONFIG.md) | The vehicle config files: fields, frame, adding a vehicle |
| [FDIR_ARCHITECTURE.md](docs/FDIR_ARCHITECTURE.md) | Modes, monitors, responses and fault injection |
| [SLOSH_MODEL.md](docs/SLOSH_MODEL.md) | The slosh model, its validation and limits |
| [KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md) | What's wrong or unfinished today, and the workarounds |
| [CONTRIBUTING.md](docs/CONTRIBUTING.md) | Branches, model rules, tests, baselines, code owners |

## Diffing and merging models

`.slx` files are binary to git, so `.gitattributes` routes them to MATLAB's merge driver (`merge=mlAutoMerge`). Set it up once per machine:
- Use your own install path (run `matlabroot` in MATLAB).
- Rerun these commands when you move to a new release.

```bash
git config --global merge.mlAutoMerge.driver '"C:/Program Files/MATLAB/<release>/bin/win64/mlAutoMerge.bat" %O %A %B %A'
git config --global difftool.mlDiff.cmd '"C:/Program Files/MATLAB/<release>/bin/win64/mlDiff.bat" $LOCAL $REMOTE'
git config --global mergetool.mlMerge.cmd '"C:/Program Files/MATLAB/<release>/bin/win64/mlMerge.bat" $BASE $LOCAL $REMOTE $MERGED'
```

- **Diffs:** `git difftool -t mlDiff -- <model>.slx` opens Simulink's Comparison Tool.
- **Conflicts:** `mlAutoMerge` resolves what it can; `git mergetool -t mlMerge` opens the three-way model merge for the rest.
- **Avoid two people editing the same `.slx`.** The shared pieces are separate subsystem files, so prefer editing one of those over the harness.

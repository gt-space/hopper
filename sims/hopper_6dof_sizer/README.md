# Hopper 6DOF Vehicle Sizer

MATLAB/Simulink sims used for preliminary vehicle sizing of the Hopper (YJSP's Propulsive Lander), modeling wind, slosh, resonant behavior, & mass/inertia variation.

**Working on the sim?** Read [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) first. It covers setup, branches, model rules, tests, and when to update the baseline.

## Quick start

Use the MATLAB release named in [`TEAM_RELEASE`](TEAM_RELEASE). CI installs that release, and models must be saved in it or older.

1. Open the project: double-click `HopperSim.prj`, or run `openProject('sims/hopper_6dof_sizer')`.
   This puts every source folder on the path and sends Simulink caches to `work/`. No `addpath` or `cd` is needed.
2. Nominal flight: `out = sim(buildSimInput())`. `buildSimInput(model, scenario)` returns a `Simulink.SimulationInput` carrying every variable the model reads, so nothing has to be in the base workspace. The old way (`sim_setup`, then simulate) is still there for the fsw_bridge scripts.
3. Other configurations: `out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('SLOSH', 'off')))` flies with slosh off. The axes: VEHICLE (hopper, ctv), CONTROLLER, SLOSH, PROPULSION, GROUND, MASS_PROPS, FDIR, FAULT_INJECTION and NAVIGATION (the nav team's GMEKF, in shadow or flown on); see [docs/VARIANTS.md](docs/VARIANTS.md).
4. Debugging: `variantConfig('SCOPES', ["attitude" "actuation"])` opens those scope groups, and `variantConfig('LOG', 'on')` writes every signal between the harness groups to a MAT file under `outputs/logs`; see [Scopes and the signal log](docs/VARIANTS.md#scopes-and-the-signal-log).
5. Monte Carlo: `R = runMonteCarlo(10)`; add `UseParallel=true` to use `parsim`. `R.summary` has one row per run. `mc_runner` and the batched/checkpointed `mc_parallelization` (for thousands of runs) use the same pieces. Every run calls CoolProp through Python (`pip install CoolProp`).

## Testing

From the project root (`buildtool` needs R2023b or later):

| Command | Runs | Time |
|---|---|---|
| `buildtool` | `check` + `test`: repo guards, unit, model-file and FDIR tests, and the nav team's filter tests | ~2 min |
| `buildtool sim` | Nominal flight against `tests/baselines/nominal.mat`, the slosh open-loop check, a smoke run of every variant choice, and the FDIR and navigation flights | ~40 min |
| `buildtool mc` | Seeded 10-run Monte Carlo against the regression limits in `tests/montecarlo/mcLimits.m` | ~5 min with a parallel pool |

### CI

[`.github/workflows/hopper-sim.yml`](../../.github/workflows/hopper-sim.yml) runs on every push and pull request that touches this folder, on GitHub-hosted Linux runners with the MATLAB release in `TEAM_RELEASE`. The repo is public, so MathWorks provides the license for these runners.

| Job | What it does |
|---|---|
| `guards` | No MATLAB, fails in seconds. `tools/ci/guards.sh`: no committed caches, no model newer than the team release. `tools/ci/codeowners_check.py`: every owner in `.github/CODEOWNERS` is a team GitHub can resolve. |
| `test` | `buildtool check test sim`. JUnit XML is uploaded as the `test-results` artifact. |
| `montecarlo` | `buildtool mc`. Writes the per-run table to the run's summary page and uploads `outputs/montecarlo/` (`mc_ci.mat`, `mc_summary.csv`). |
| `model-diff` | Pull requests only. A Simulink Comparison HTML report for every changed `.slx`, uploaded as the `model-diffs` artifact. |

To make `main` require these jobs, an admin turns this on once under **Settings → Branches → Add branch protection rule** for `main`:
- **Require status checks to pass:** `Repo guards`, `Unit, model and nominal-flight tests`, `Monte Carlo (10 seeded runs)`.
- **Require a pull request before merging,** with 1 approval.

What each task checks:

- **`check`** fails if build caches or autosaves are committed, or if a model is saved in a release newer than `TEAM_RELEASE`.
- **`mc`** fails a run that errors, diverges, or does worse than today's vehicle plus margin. Separately, `R.summary.meets_requirements` reports each run against the stricter mission requirements. Edit those in the MISSION CONSTRAINTS section of `inputs/mission_inputs.m` (fields marked `[req]`).
- **Intended physics change:** if your change is meant to alter the flight, run `NominalFlightTest.updateBaseline()`, say so in the PR, and commit the new baseline.

Results are repeatable: the same inputs give bit-identical results in any session.

Save models in the **team release** (`TEAM_RELEASE`) or older. A model saved in a newer release can't be opened by everyone else.

## Layout

| Folder | Contents |
|---|---|
| `models/` | `hopper_6dof.slx` is the master model, a variant harness (see `docs/VARIANTS.md`); `hopper_6dof_NED_v2_fswBridge.slx` is the flight-software bridge. `dev/` holds work in progress (LQI); `standalone/` holds component models. v2 and its feature copies are in `archive/models/`. |
| `subsystems/` | Shared subsystem files the models reference: slosh (`Slosh` and the `Slosh_Axis` pendulum used six times), mass_properties, actuators (TVC/RCS, EDF fans), environment (wind), ground (#27's strut ground contact), gnc (LQR gain scheduler), propulsion (mdot lookup, propellant mass, sensors, throttle valve controller), flight_logic (engine cutoffs, touchdown stop), FDIR, navigation filters. Edit a block here once and every model using it changes; see CONTRIBUTING, *Shared subsystems*. |
| `src/` | MATLAB code: `setup`, `sizing` (design trades with `sizeVehicle`; not run in sim setup), `dynamics`, `gnc`, `propulsion`, `slosh`, `sensors` (the shared sensor model FDIR and navigation read), `navigation` (the adapter between the sim and the nav team's filter), `fdir`, `faults` (fault injection), `monte_carlo`, `analysis`, `util`. |
| `vehicles/` | One config file per vehicle (`hopper/hopper_config.m`, `ctv/ctv_config.m`): geometry, masses and hardware, the sim's source of truth for the vehicle. See `docs/VEHICLE_CONFIG.md`. |
| `inputs/` | `mission_inputs.m` (requirements, trajectory, environment, dispersions), `ctv_inputs.m` and `mc_params.json`. |
| `data/` | Lookup tables the sims read (`mdot_lookup.xlsx`, `wind_vectors2.mat`, `cg_I_LUT.mat`, `baffle_lut.*`, `engine_contour.xlsx`). |
| `results/` | Tracked Monte Carlo and sweep results (`.mat`). |
| `tests/`, `validation/` | Automated tests and validation studies. |
| `fsw_bridge/` | C codegen and Rust SIL bridge to flight software. |
| `tools/refactor/` | Scripts that change the harness or shared subsystem files reproducibly (for example `parameterizeVehicleGeometry`). |
| `docs/` | Model notes: `VEHICLE_CONFIG.md`, `VARIANTS.md`, `FDIR_ARCHITECTURE.md`, `SLOSH_MODEL.md`. |
| `archive/` | Superseded models and scripts, kept for reference and **not** on the path. |
| `outputs/`, `work/` | Generated figures, videos, CSVs and build caches. Gitignored. |

If you move or add files with `git mv` or Explorer instead of the Project panel, run `syncProjectFiles` so the project's file list matches the folders.

Scripts that write files use `hopperRoot('data', ...)`, `hopperRoot('results', ...)` or `hopperRoot('outputs', ...)`, so output lands in the right folder whatever the current folder is.

## Diffing and merging models

`.slx` files are binary to git. `.gitattributes` routes them to MATLAB's merge driver (`merge=mlAutoMerge`). Set it up once per machine. Use your own install path (run `matlabroot` in MATLAB), and rerun these when you move to a new release:

```bash
git config --global merge.mlAutoMerge.driver '"C:/Program Files/MATLAB/<release>/bin/win64/mlAutoMerge.bat" %O %A %B %A'
git config --global difftool.mlDiff.cmd '"C:/Program Files/MATLAB/<release>/bin/win64/mlDiff.bat" $LOCAL $REMOTE'
git config --global mergetool.mlMerge.cmd '"C:/Program Files/MATLAB/<release>/bin/win64/mlMerge.bat" $BASE $LOCAL $REMOTE $MERGED'
```

- **Diffs:** `git difftool -t mlDiff -- <model>.slx` opens Simulink's Comparison Tool.
- **Conflicts:** `mlAutoMerge` resolves what it can. For the rest, `git mergetool -t mlMerge` opens the three-way model merge.

Avoid having two people edit the same `.slx`. The shared pieces live in separate subsystem files, so changes rarely overlap; prefer editing one of those over the top-level model.

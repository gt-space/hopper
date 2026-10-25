# Hopper 6DOF Vehicle Sizer

MATLAB/Simulink sims used for preliminary vehicle sizing of the Hopper (YJSP's Propulsive Lander), modeling wind, slosh, resonant behavior, & mass/inertia variation.

**Working on the sim?** Read [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) first. It covers setup, branches, model rules, tests, and when to update the baseline.

## Quick start

Use the MATLAB release named in [`TEAM_RELEASE`](TEAM_RELEASE). CI installs that release, and models must be saved in it or older.

1. Open the project: double-click `HopperSim.prj`, or run `openProject('sims/hopper_6dof_sizer')`.
   This puts every source folder on the path and sends Simulink caches to `work/`. No `addpath` or `cd` is needed.
2. Nominal flight: `out = sim(buildSimInput())`. `buildSimInput(model, scenario)` returns a `Simulink.SimulationInput` carrying every variable the model reads, so nothing has to be in the base workspace. The old way (`sim_setup`, then simulate `hopper_6dof_NED_v2`) still works.
3. Monte Carlo: `R = runMonteCarlo(10)`; add `UseParallel=true` to use `parsim`. `R.summary` has one row per run. `mc_runner` and the batched/checkpointed `mc_parallelization` (for thousands of runs) use the same pieces. Every run calls CoolProp through Python (`pip install CoolProp`).

## Testing

From the project root (`buildtool` needs R2023b or later):

| Command | Runs | Time |
|---|---|---|
| `buildtool` | `check` + `test`: repo guards, unit, model-file and FDIR tests | ~2 min |
| `buildtool sim` | Nominal flight against `tests/baselines/nominal_v2.mat`, plus the slosh open-loop check | ~5 min |
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
| `models/` | Top-level models (`hopper_6dof_NED_v2.slx` is the master). `dev/` holds work in progress (LQI); `standalone/` holds component models. |
| `subsystems/` | Subsystem files referenced by the models (FDIR, actuators, propulsion sensors/valves, navigation filters). |
| `src/` | MATLAB code: `setup`, `sizing`, `dynamics`, `gnc`, `propulsion`, `slosh`, `fdir`, `monte_carlo`, `analysis`, `util`. |
| `inputs/` | `mission_inputs.m` and `mc_params.json`. |
| `data/` | Lookup tables the sims read (`mdot_lookup.xlsx`, `wind_vectors2.mat`, `cg_I_LUT.mat`, `baffle_lut.*`, `engine_contour.xlsx`). |
| `results/` | Tracked Monte Carlo and sweep results (`.mat`). |
| `tests/`, `validation/` | Automated tests and validation studies. |
| `fsw_bridge/` | C codegen and Rust SIL bridge to flight software. |
| `tools/refactor/` | Scripts that build or convert models reproducibly (for example `build_v2_fluid_integration`). |
| `docs/` | Model notes (`SLOSH_MODEL.md`). |
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

Avoid having two people edit the same `.slx`. Shared pieces are moving into separate subsystem files so changes rarely overlap.

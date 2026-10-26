# Contributing to the hopper 6DOF sim

How we work on this sim together without breaking each other's models. The [README](../README.md) covers the folder layout and how to run things.

## One-time setup

1. **MATLAB, in the release named in [`TEAM_RELEASE`](../TEAM_RELEASE).** Toolboxes:
   - Simulink, Aerospace Blockset
   - Communications, DSP System and Signal Processing Toolbox
   - Simscape, Simscape Fluids
   - Control System, Optimization, and Statistics and Machine Learning Toolbox
   - Parallel Computing Toolbox (optional, for parallel Monte Carlo)

   Other releases can be installed alongside, but only **save** from the team release (see [Models](#models)).
2. **Python with CoolProp** (`pip install CoolProp`). Every Monte Carlo run calls it. Check with `pyenv` in MATLAB.
3. **Git merge tools for `.slx`:** run the three `git config` lines in the README's *Diffing and merging models* section.
4. Open `HopperSim.prj` every time you work. It sets the path; there is no `addpath`.

## Day-to-day workflow

- **Branch from `main`, keep it short.** One feature or fix per branch, merged within a week or two.
  - Merge `main` into your branch at least weekly.
  - Branches that drift for a month end up with `.slx` conflicts nobody can resolve.
- **Say which model you are editing** in the team channel before changing a shared `.slx`. Two people editing the same binary file at once means one of them redoes their work.
- **Never copy the top-level model to add a feature.** Every `_ground`, `_fluid_integration` or `_fault_injection` copy so far drifted from the master within weeks.
  - Put the new piece in its own subsystem file under `subsystems/<area>/`.
  - Wire it into the master model behind a switch.
  - Features that replace a piece of the vehicle become a choice in the variant harness `models/hopper_6dof.slx`; see [VARIANTS.md](VARIANTS.md). Ask if unsure where yours goes.
- **Open a pull request early** (draft is fine). CI runs on every push, and the PR template asks for the Monte Carlo numbers.
- **Code owners** (`.github/CODEOWNERS`) are GitHub teams. A pull request automatically requests review from the team that owns each file it touches:
  - `hopper-gca`: guidance and control
  - `hopper-nav`: navigation filters
  - `hopper-fdir`: FDIR
  - `hopper-fsw`: FSW bridge
  - `hopper-modsim`: everything else, including propulsion and structures for now
  - `hopper-gnc-leads`: the files that set the rules for everyone. It is the sole owner of `TEAM_RELEASE`, `tests/montecarlo/mcLimits.m`, the CI workflow and CODEOWNERS itself, and co-owns `inputs/mission_inputs.m`, where the requirements live.

  Joining or leaving only changes team membership (gt-space org → Teams), never the file. Edit CODEOWNERS only when an area moves to a different team or a new team is created. CI's guards job fails if GitHub can't resolve an owner, for example a renamed team or one without Write access to the repo.

## Files and folders

- **Adding or moving files:**
  - Easiest: the MATLAB Project panel.
  - With `git mv` or Explorer: run `syncProjectFiles` afterwards.
  - New folders need `addPath(currentProject, folder)` once.
- **Scripts that write files** use `hopperRoot`:
  - lookup tables → `hopperRoot('data', ...)`
  - tracked results → `hopperRoot('results', ...)`
  - figures, videos, CSVs → `hopperRoot('outputs', ...)`, which is gitignored

  Never write to the current folder.
- **Superseded code** goes to `archive/`, not the trash. `archive/` is off the path, so it can't shadow live code. Delete only obvious junk.
- **Don't commit generated files.** `slprj/`, `*.slxc`, `*.asv`, `.slx.r20xx` backups and anything in `outputs/` or `work/` are ignored, and CI's guards job fails if one sneaks in.

## Models

- **Save in the team release or older.** A model saved in a newer release can't be opened by anyone on the team release, and CI rejects it.
  - If you saved one by accident, re-save it from the newer release with `Simulink.exportToVersion(model, 'tmp.slx', teamRelease())` and replace the file.
- **Give blocks real names:** `Thrust_to_Mdot_LUT`, not `Subsystem3`. Scripts and tests find blocks by path, so renaming later breaks them.
- **Read parameters from workspace variables or `IN.*`, not hard-coded literals.** Then `buildSimInput` can set them per run.
- **Random blocks need a fixed seed.** For example, Bernoulli/Random Source `SeedSource = 'Parameter'`. A block on `'Auto'` reads MATLAB's global random stream and makes every run depend on what ran before it.
- **Reviewing someone's model change:** download the PR's **model-diffs** artifact (a Simulink Comparison report per changed `.slx`), or run `git difftool -t mlDiff -- <file>.slx` locally.

## Shared subsystems

Blocks that more than one model uses live in their own files under `subsystems/<area>/`, for example `subsystems/propulsion/Thrust_to_Mdot_LUT.slx`. The models hold **Subsystem Reference** blocks that point at those files. The reference block shows the file name under its icon.

- **To change a shared block, open its file:** double-click the reference block, then edit and save. Every model that references the file picks up the change. Prefer this to editing a copy inside one model.
- **One file can serve several blocks.** `Slosh_Axis.slx` is all six slosh pendulums (ox/fuel × axial/pitch/yaw). Each instance gets its tank and axis only from its wiring, so keep tank- or axis-specific names out of the file.
- **Settings that differ per use are mask parameters.** The fuel and ox valves share `propulsion/valves/Throttle_Valve_Controller.slx`; each instance sets its own `ramp_slope` in the block's mask dialog. Add a mask parameter rather than a second copy of the file.
- **Files can nest.** `Slosh.slx` references `Slosh_Axis.slx`, and `Actuators_TVC_RCS.slx` references `EDF_Dynamics.slx`.
- **Tests that tweak a shared block work on a copy.** `set_param` on a block inside a reference edits the shared file. Copy the file's contents into a plain subsystem first, as `validation/slosh_openloop/sloshOpenLoopTest.m` does with `Simulink.BlockDiagram.copyContentsToSubsystem`.
- **Keep the ports stable.** Port names and order are the interface every referencing model is wired to. If you must add or reorder a port, update every model that uses the file in the same PR. CI's model tests load them all.
- **Read parameters from workspace variables** (for example `IN.*` or `thrust_data`), not from the parent model. The file has no workspace of its own, and `buildSimInput` sets the variables per run.
- **To share another block,** add it to the list in `tools/refactor/convertToSubsystemReferences.m` and run it. It converts the block in v2, and models built from v2 (such as the fluid-integration model) pick it up when rebuilt. Check that the nominal flight is unchanged: `buildtool sim`.
- **Conflicts get rarer.** Two people editing different subsystem files never touch the same `.slx`, so prefer a new or existing subsystem file over editing the top-level model.

## Moving to a new MATLAB release

The team release lives in one place: [`TEAM_RELEASE`](../TEAM_RELEASE), one line such as `R2026a`. CI installs that release, and `tools/ci/guards.sh`, `buildtool check` and `tests/models` all check models against it (`teamRelease()` in MATLAB). Models saved in older releases stay valid, so upgrading never forces anyone to re-save files.

To move up, in one PR of its own:

1. **Agree on a date with the team.** Everyone needs the new release installed before the PR merges, or they can't open models saved after it. `hopper-gnc-leads` owns `TEAM_RELEASE` and approves the PR.
2. **Change the one line in `TEAM_RELEASE`.**
3. **Run `buildtool check test sim mc` locally on the new release.**
   - A new solver or block version can shift the nominal flight by tiny amounts. If `NominalFlightTest` fails only by small numbers, run `NominalFlightTest.updateBaseline()` in the same PR and say so.
   - A large change means the release changed behavior: find out why before merging.
4. **Check the CI run.** It installs the new release on its own.
   - If `setup-matlab` rejects a product name, fix `MATLAB_PRODUCTS` in `.github/workflows/hopper-sim.yml`.
   - If MATLAB rejects the Python version, change `python-version` in that file.
5. **After merging,** everyone reruns the three `git config` lines in the README with the new install path, so model diffs and merges use the new release.

## Running and testing

| What | How |
|---|---|
| One nominal flight | `out = sim(buildSimInput())` |
| One scenario | `scenarios = generateScenarios('mc_params.json', 10); out = sim(buildSimInput('hopper_6dof_NED_v2', scenarios(3)))` |
| Monte Carlo | `R = runMonteCarlo(100, UseParallel=true)`, then look at `R.summary` |
| Before pushing | `buildtool` (about 2 min), plus `buildtool sim` if you touched the model or setup |

`runMonteCarlo` sizes its parallel pool from free RAM, since each worker takes about 3.3 GB. Pass `Workers=N` to override it, but don't go above what fits: Windows will page to disk, and the page file can grow by tens of GB.

### When the flight is supposed to change

`buildtool sim` compares the nominal flight with `tests/baselines/nominal_v2.mat`. If your change is meant to alter the flight (new gains, new physics):

1. Run `NominalFlightTest.updateBaseline()` and commit the new `.mat`.
2. Say so in the PR, with the Monte Carlo before/after numbers.

If you did **not** mean to change the flight and `sim` fails, that's a regression: find it, don't re-baseline it.

### Requirements and CI limits

- **Mission requirements** live in the *MISSION CONSTRAINTS* section of `inputs/mission_inputs.m` (fields marked `[req]`). Each Monte Carlo run reports `meets_requirements`.
- **CI regression limits** live in `tests/montecarlo/mcLimits.m`. They describe today's vehicle plus margin, and the CI Monte Carlo fails outside them.

Tighten the regression limits as the vehicle improves. Loosening one needs a reason in the PR and an approval from `hopper-gnc-leads`, which owns `mcLimits.m`.

### Adding tests

`matlab.unittest` classes:
- `tests/unit/`: fast, no simulation
- `tests/models/`: model files load and are saved in the team release or older
- `tests/sim/`: a single simulation
- `tests/montecarlo/`: the Monte Carlo

`buildfile.m` picks up new test files in those folders automatically.

## Commit messages

First line: what changed, in the imperative ("Add strut damping to ground contact"). Body: why, and what you checked. For a model change, name the subsystems touched; the binary diff can't say it for you.

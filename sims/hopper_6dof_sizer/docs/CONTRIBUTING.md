# Contributing to the hopper 6DOF sim

How we work on this sim together without breaking each other's models. The [README](../README.md) covers the folder layout and how to run things.

## One-time setup

1. **MATLAB R2026a.** Toolboxes:
   - Simulink, Aerospace Blockset
   - Communications, DSP System and Signal Processing Toolbox
   - Simscape, Simscape Fluids
   - Control System, Optimization, and Statistics and Machine Learning Toolbox
   - Parallel Computing Toolbox (optional, for parallel Monte Carlo)

   Other releases can be installed alongside, but only **save** from R2026a (see [Models](#models)).
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
  - The shared-subsystem/variant work is turning these switches into proper variant choices; ask if unsure where yours goes.
- **Open a pull request early** (draft is fine). CI runs on every push, and the PR template asks for the Monte Carlo numbers.
- **Code owners** (`.github/CODEOWNERS`) are asked to review changes in their area. If you now own an area, add yourself there.

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

- **Save in R2026a or older.** A model saved in a newer release can't be opened by anyone on R2026a, and CI rejects it.
  - If you saved one by accident, re-save it from the newer release with `Simulink.exportToVersion(model, 'tmp.slx', 'R2026a')` and replace the file.
- **Give blocks real names:** `Thrust_to_Mdot_LUT`, not `Subsystem3`. Scripts and tests find blocks by path, so renaming later breaks them.
- **Read parameters from workspace variables or `IN.*`, not hard-coded literals.** Then `buildSimInput` can set them per run.
- **Random blocks need a fixed seed.** For example, Bernoulli/Random Source `SeedSource = 'Parameter'`. A block on `'Auto'` reads MATLAB's global random stream and makes every run depend on what ran before it.
- **Reviewing someone's model change:** download the PR's **model-diffs** artifact (a Simulink Comparison report per changed `.slx`), or run `git difftool -t mlDiff -- <file>.slx` locally.

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

Tighten the regression limits as the vehicle improves. Loosening one needs a reason in the PR.

### Adding tests

`matlab.unittest` classes:
- `tests/unit/`: fast, no simulation
- `tests/models/`: model files load and are saved in R2026a or older
- `tests/sim/`: a single simulation
- `tests/montecarlo/`: the Monte Carlo

`buildfile.m` picks up new test files in those folders automatically.

## Commit messages

First line: what changed, in the imperative ("Add strut damping to ground contact"). Body: why, and what you checked. For a model change, name the subsystems touched; the binary diff can't say it for you.

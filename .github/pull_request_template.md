<!-- Hopper 6DOF sim? See sims/hopper_6dof_sizer/docs/CONTRIBUTING.md. Delete sections that don't apply to your PR. -->

## What and why

<!-- One or two sentences: what this changes and why. Link the issue or Notion page if there is one. -->

## Sim changes

- [ ] No `.slx` changed
- [ ] Changed models/subsystems (list them; CI attaches a **model-diffs** report for each):
  -

Saved in the team release (`sims/hopper_6dof_sizer/TEAM_RELEASE`) or older? <!-- CI's guards job fails otherwise -->

Copying the top-level model to add a feature makes merges painful. If you made a new variant this way, say so here so it can become a variant choice instead.

Changed a file in `subsystems/navigation` (the nav team's)? Say which, and tell them: they merge from their own branch.

## Configurations

Vehicle and variant choices this affects or was run with (`variantConfig`, for example `VEHICLE=ctv` or `NAVIGATION=gmekf`; "defaults" for the nominal hopper):
-

## Results

- [ ] **No intended change to the flight**: the nominal baseline and the Monte Carlo summary should match `main`
- [ ] **Intended physics/controls change**: I ran `NominalFlightTest.updateBaseline()` and committed the new `tests/baselines/nominal.mat`

Monte Carlo (from the CI run's summary page or `buildtool mc`: 10 hopper runs, 5 CTV). Paste what moved:

| | hopper `main` | hopper this PR | CTV `main` | CTV this PR |
|---|---|---|---|---|
| Within CI limits | /10 | /10 | /5 | /5 |
| Meets mission requirements | /10 | /10 | /5 | /5 |
| Impact speed (max, m/s) | | | | |
| Landing offset (max, m) | | | | |
| Max tilt (deg) | | | | |

Touched FDIR, navigation or the controller? Say which of these flights you checked and what moved (`buildtool sim` runs them all):
- FDIR fault cases (`tests/sim/FdirFlightTest.m`)
- Navigation in shadow or flown on the estimate, hopper and CTV (`tests/sim/NavigationFlightTest.m`)

Changed limits or requirements? <!-- tests/montecarlo/mcLimits.m or the [req] fields in inputs/mission_inputs.m; say why -->

## Tested

- [ ] `buildtool` (check + test, ~2 min) passes locally
- [ ] `buildtool sim` passes locally, if you changed a model, a subsystem or setup code (~40 min; CI runs it on every push either way)
- Other variants or models run by hand:

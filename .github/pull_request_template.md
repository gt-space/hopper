<!-- Hopper 6DOF sim? See sims/hopper_6dof_sizer/docs/CONTRIBUTING.md. Delete sections that don't apply to your PR. -->

## What and why

<!-- One or two sentences: what this changes and why. Link the issue or Notion page if there is one. -->

## Sim changes

- [ ] No `.slx` changed
- [ ] Changed models/subsystems (list them; CI attaches a **model-diffs** report for each):
  -

Saved in **R2026a** or older? <!-- CI's guards job fails otherwise -->

Copying the top-level model to add a feature makes merges painful. If you made a new variant this way, say so here so it can become a variant choice instead.

## Results

- [ ] **No intended change to the flight**: the nominal baseline and the Monte Carlo summary should match `main`
- [ ] **Intended physics/controls change**: I ran `NominalFlightTest.updateBaseline()` and committed the new `tests/baselines/nominal_v2.mat`

Monte Carlo (from the CI run's summary page or `buildtool mc`). Paste what moved:

| | `main` | this PR |
|---|---|---|
| Within CI limits | /10 | /10 |
| Meets mission requirements | /10 | /10 |
| Impact speed (max, m/s) | | |
| Landing offset (max, m) | | |
| Max tilt (deg) | | |

Changed limits or requirements? <!-- tests/montecarlo/mcLimits.m or the [req] fields in inputs/mission_inputs.m; say why -->

## Tested

- [ ] `buildtool` (check + test) passes locally
- [ ] `buildtool sim` passes locally
- Other variants or models run by hand:

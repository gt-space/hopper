# Onboard G-FOLD (archived)

Not used by the CTV sim. Kept so we can pick it back up next semester.

The team decided not to fly G-FOLD guidance on CTV for now. The CTV sim flies
the precomputed reference trajectory with LQR gains computed offline
(`ctv/LinerizationMaster_ctv.m`), and the controller runs those gains in flight.

## What's here

- `gfold_onboard.m`: fuel-optimal descent solve from the current state (CVX, or `coneprog` as a fallback). Limits come from `mission_inputs()`.
- `gfold_onboard_step.m`: stateful re-planner. Re-solves at most once per second and splices the result into the reference table. Raises flag F-23 if the solve fails.
- `gfold_replan.m`: body of the Simulink MATLAB Function block.
- `add_onboard_gfold.m`: inserts the `GFOLD Replan` block into the model through the Simulink API. **This edits the .slx.** Don't run it unless you're deliberately picking this work back up.

None of this has been run in MATLAB yet.

The offline planners in `sims/GFOLD/` (`gfold_cvx.m`, `gfold_descent_v3.m`) are separate and stay where they are.

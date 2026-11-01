# Known issues

What's wrong, unfinished or surprising in the sim today (2026-10-31), and what to do about it. Fix one? Update this page in the same pull request.

## Vehicle performance

- **Most dispersed hopper runs miss the landing-offset requirement.**
  - In CI's 10 seeded runs (`buildtool mc`), all 10 pass the regression limits, but only 1 meets every mission requirement.
  - The other 9 land 1.3–2.1 m from the pad against the 1 m requirement (`IN.mission.max_landing_offset`). Every other requirement holds.
  - Re-fly one case with `R.scenarios(k)` and `SCOPES` to investigate ([User guide](USER_GUIDE.md#debugging-a-monte-carlo-run)).
- **The nominal hopper lands 0.66 m off the pad.** That's inside the 1 m requirement, but it was 0.43 m on the earlier stand-in mass properties. The CAD dry mass properties (114.894 kg, since 2026-10-29) moved it.
- **The dispersed CTV mostly fails.**
  - The nominal CTV meets every requirement in calm air and lands on the nav estimate.
  - In `CtvMonteCarloTest`'s 5 seeded runs, only 1 meets requirements; the rest miss the landing offset, one also the touchdown tilt. Across 10 seeded runs, most land metres off or tilt badly, and some diverge.
  - Suspects, in order: LQR tuning on hopper-derived weights, products-of-inertia coupling (dropped from the flown tensor), and the unmodelled EDF reaction torque. See [VARIANTS: The CTV](VARIANTS.md#the-ctv).
- **`CONTROLLER=lqi` is experimental.** LQI-DEV's fixed-gain hover regulator was designed on a toy plant and isn't tuned. It compiles and runs, and the variant test only checks that. It may not hold the vehicle up.

## Model fidelity

- **Structural bending is a stopgap.** The hopper's first bending mode couples into the TVC through `k_engine = (engine cg_x / height)^2` with a 10 Hz, 1 % damped mode (`VC.structure.modal`). The numbers are placeholders for a real structural model. The CTV has it switched off (`modal.enabled = false`).
- **Sizer numbers still in the hopper config.** The tank wall masses and leg stiffness and mass are the conceptual sizer's 2026-10 values (tagged `sizer 2026-10` in `vehicles/hopper/hopper_config.m`), not as-built. `importCadMassProps` is a stub; enter CAD values by hand for now.
- **The CTV has many unconfirmed values.** Values marked `ASSUMPTION` in `vehicles/ctv/ctv_config.m` and `inputs/ctv_inputs.m`:
  - EDF and gimbal heights and the EDF mass
  - throttle range
  - roll-prop thrust and arm (still the hopper's tank radius)
  - leg spread, aero, ascent and descent limits

  Also missing from the CTV: the lateral CG (kept in the config but not flown), the products of inertia, and the EDF reaction torque.
- **Slosh assumes flat-bottom tanks and a single mode per tank.** Domed bottoms would change exactly the low-fill behaviour at landing. See [SLOSH_MODEL.md §9](SLOSH_MODEL.md#9-known-limitations).
- **Sensor noise isn't from data sheets yet.** The noise comes from `sim_sensor_params`, the values FDIR's thresholds were tuned on, not the ADIS16500 data sheet. The lidar mount and GPS antenna positions are the nav team's placeholders (`ASSUMPTION` in the config).

## FDIR and navigation

- **FDIR thresholds and persistence times are placeholders.** They are to be set from Monte Carlo before PDR.
- **The engine cutoff differs from the PDR.** The sim cuts thrust at 0.05 m (`IN.mission.engine_cutoff_alt`); the PDR says 0.1 m.
- **The FMECA rows outside GNC aren't written.**
- **A stuck gyro with `NAVIGATION=gmekf` still crashes.** FDIR masks the gyro and the filter stops fusing it. Without a gyro the GMEKF can't hold attitude on differential GPS and lidar alone, and the abort flies on that attitude.
- **The nav team's newest commits aren't merged.** `c32fba9`..`69703bc` ("flying the gyro", work in progress) are on their branch only. `tools/refactor/syncNavigation.m` merges them when they're ready. Several of their files carry changes made here; the nav team note (`outputs/nav_team_note.md`, local) lists them.

## Tooling and workflow

- **`y = logical ...` prints to the console** many times per run. It's harmless, and the source hasn't been found yet. It's most likely an expression missing its semicolon in code that runs while the model compiles or simulates; bisecting by variant choice is the next step.
- **The Run button doesn't work on a freshly opened harness.** The model's variables come from `buildSimInput`, not the base workspace.
  - Workaround: `buildSimInput(...).applyToModel()`, then press Run. Close the model without saving afterwards ([User guide](USER_GUIDE.md#running-from-the-simulink-window)).
  - `sim_setup` predates the variants and doesn't set their controls.
- **A newer MATLAB on the PATH silently saves models in the wrong release.** On machines with several releases, `matlab` from a terminal may be R2026b. A model saved from it fails CI's release check and leaves a `.slx.r2026a` backup. Start R2026a explicitly for anything that saves.
- **`throttle_valve_simulink.slx` was saved in R2026b.** R2026a exports it on load, which prints a notice in the tests. Its owners need to re-save it from R2026a (`Simulink.exportToVersion`).
- **The top-level layout is automatic.** `arrangeSystem` places the seven groups on a diagonal. No line crosses a block, but it isn't hand-arranged.
- **The signal log is large.** About 42 MB per hopper flight, mostly `nav_bus.p_cov`, with no decimation option. Keep `LOG` off in big Monte Carlos.
- **Scope groups haven't been exercised in the GUI.** `VariantHarnessTest` builds and flies all nine with the log on, but nobody has yet checked how their windows look interactively.
- **`buildtool sim` takes ~40 minutes,** because it flies every variant choice and the long FDIR and navigation cases. Run `buildtool check test` (~2 min) on every change, and the sim tests when the model or GNC changes.
- **Monte Carlo needs CoolProp.** `prop_system` calls Python CoolProp in every run. A single flight doesn't need it.

## Refactor status

- **Not merged yet.** Everything here is on `refactor/block-names` (worktree `hopper-integrate`), unpushed. Merging is pending the team's review.
- **Repo settings still to do.** An admin needs to set up branch protection (steps in the [README](../README.md#ci)) and the GitHub teams named in CODEOWNERS ([CONTRIBUTING](CONTRIBUTING.md#day-to-day-workflow)).

# Hopper FullFlow design and progress

This is the working design document for the propulsion simulator. Update the
checklist and validation notes when behavior changes. Inputs belong in
`HOPPER_DATA.py`; generated results belong in `Hopper_Results/`.

## Objective

Given a thrust history and configurable hardware, predict feed-system pressures,
temperatures, flow, inventories, and engine performance over time.

## Current implementation work

- [x] Pressure versus component station, at selectable simulation times.
- [x] Prescribed fuel-regen temperature rise versus chamber pressure.
- [x] Automatically generated component/source box diagram.
- [x] Selectable flight or parameterized ramp profile; demonstrated maximum report.
- [x] Separate valve opening/closing travel times and controller sampling.
- [x] Propellant residual masses, consumed masses, and percentages.
- [x] Explicit tank pressure deadband inputs.
- [x] Individual COPV states and minimum-temperature flags.
- [x] Engine efficiencies from functions or CSV tables.
- [x] Automated checks and integration/plot validation.

## Agreed definitions

- A COPV minimum temperature is an operating-limit check. Crossing it flags the
  run; it does not clamp temperature or stop flow.
- Maximum thrust means the highest sustained thrust demonstrated during the
  selected run. A run that never reaches a capacity limit provides a lower bound,
  not a proven global maximum. Report this distinction explicitly.
- Residual means remaining liquid in each tank. Feed-line trapped inventory is
  outside the current zero-volume liquid-junction model.
- Regen is prescribed temperature rise, not a wall/combustion heat-transfer model.
  Its effect on downstream liquid properties is included; wall thermal inertia
  and feedback into equilibrium combustion temperature are not implied.
- Valve travel times mean full-stroke times. Opening and closing may differ.
  Partial strokes scale with travel distance. Zero means ideal instantaneous.
- Controller sample period is independent of valve travel time. Solver resolution
  must still be fine enough to resolve both.

## Validation

Twelve focused behavioral tests pass. All 19 named plots rendered with CSV output;
the component diagram and station snapshots were visually checked. A complete
1.5 s startup/burn/shutdown check passed at 5 ms integration/output intervals,
using a custom 500-to-700 N demand, nonzero regen rise, and an intentionally
violated COPV temperature limit. This was a feature check, not a full flight
profile rerun or hardware calibration.

`tests/check_integration.py` checks a saved feature-validation run for conservation,
regen temperatures, valve travel rates, limits, report behavior, and property
interpolation accuracy. Tank inventory matches integrated tank outflow. Total
pressurant mass drift was below 0.01 g in the checked run. At discrete events,
saved connector rates are endpoint samples; inventory conservation is checked
over the entire connected gas system. Checked property interpolation error was
below 0.003% relative to direct ThermoProp queries.

Repeat with `python -m unittest discover -s tests -v` and
`python tests/run_feature_validation.py`. The latter applies temporary in-process
overrides and writes to `validation_outputs/`. Existing user-generated results
are preserved; new diagram/station files were added to `Hopper_Results/` from the
current configuration and existing history respectively.

## Implementation decisions

- `Hopper_Inputs.py`: input validation, ramp sampling, function/CSV efficiencies,
  and interpolation utilities. CSV format is in the README and `inputs/` examples.
- `Hopper_Control.py`: asymmetric linear valve travel based only on accepted
  state history (safe across retries), and sampled pressure hysteresis. Sample
  instants are solver breakpoints. Physical response is distinct from controller
  sampling and command-demand smoothing.
  Its `FunctionLookup` adapter disables FullFlow's object reuse for functions
  returning dictionaries. Otherwise `dict.update` can leave calculated outputs
  frozen at their initial values. Pressure-dependent efficiency behavior has a
  regression test, and regen behavior is checked against the coupled run.
- `Regen` in `Hopper_Elements.py`: hydraulic restriction plus prescribed dT(Pc).
  Downstream fluid properties use a precomputed temperature table to avoid
  expensive property calls inside nonlinear iterations. Increase
  `REGEN_PROPERTY_TEMPERATURE_POINTS` and check convergence if needed.
- COPVs have independent mass/energy states and configurable outlet restrictions
  into a common header with positive configurable gas volume. A zero-volume header
  was underdetermined when all delivery valves closed; retaining trapped gas
  gives it a physical pressure and energy state. Its 0.05 L default is a placeholder.
  Legacy combined pressure/temperature tracks now describe
  the header; combined mass is the bottle total.
- `Hopper_Report.py`: records run-specific timing, inventories, throttle limits,
  and COPV thresholds in HDF5. Reports use those saved inputs. Old histories lack
  the metadata needed for the new summary and must be rerun.
- `Hopper_Diagram.py`: derives a PNG and scalable SVG from configured components.
  Station plots identify the source and each component outlet. Gas-group pressure
  allocation remains an approximation and is labeled accordingly.
- Default efficiencies now use bounded LOW-to-HIGH interpolation. The former
  `max(1, ...)` expression did not implement the intended interpolation. Existing
  comparisons are historical and must not be treated as current calibration.

## Model boundaries

The chamber remains quasi-steady and zero-volume. Liquid feeds remain single phase.
Measured geometry, valve timing, regen curves, efficiencies, and COPV temperature
limits must be supplied by the user; example/default values are not hardware data.

# Hopper FullFlow design and progress

This is the working design document for the propulsion simulator. Update the
checklist and validation notes when behavior changes. Inputs belong in
`hopper_config.json`; generated results belong in timestamped folders under
`Hopper_Results/`. `HOPPER_DATA.py` is now a compatibility adapter, not an input file.
See [HANDOFF.md](HANDOFF.md) for the complete delivered-commit history, final
checked-in inputs, and consolidated validation/limitations. Earlier sections
below retain development-stage input and test-count history.

## November 7, 2026: merging the element-list model into the JSON version

The JSON replacement (`996f8509`) was built from the October 26 commit, so it left
out the work from October 27 to November 3 (`63ef5f4`). That work is now ported into
the JSON structure. Inputs: where only one side changed a value since October 26 that
value stands; where both changed it, the element-list value stands. Notes in the
config record each source.

Brought back from `63ef5f4`: the `Backend/`, `isolated_sims/`, `thrust_curves/`,
`Outputs/` and `Old_sims_comparison/` layout; collapse factor in place of ullage HTC;
side paths, regulator, reliefs, purges, drain and torch; thrust-curve CSV loading; the
180 lbf/s slew limit and ramp breakpoints; efficiencies sampled at the previous step
in the transient; the run report, component table, schematic and parameter workbook;
vent, drain, tank-sweep and purge calculators; the Monte Carlo; `ASSUMPTIONS.txt`.

Kept from the JSON version: the JSON loader and validation, run archives, `Regen`
heating, selectable efficiency modes, station, regen and COPV-temperature plots, the
component diagram, sustained-thrust and residual reports, valve travel, sampled
bang-bang, the bounded nonlinear solve, and the tests.

Decisions:

- COPVs are lumped by default (`pressurization.copv_model`). The bottle-and-header
  model cannot converge at the 0.05 s step: its first step stalls. Identical bottles
  give the same answer lumped. Header mode stays available for bottles that differ.
- Valve travel times default to 0 and the controller samples at 0.05 s, which is the
  element-list behaviour. Nonzero travel caps the step at 1/10 of a stroke.
- `HOPPER_DATA.py` names every input statically instead of `globals().update(...)`, so
  editors resolve `from HOPPER_DATA import *`. It raises if a config value is unnamed.
- JSON notes replace the `#EDIT` comments the table and workbook used to read.
- Fixed a header energy bug: `Volume.total_enthalpy_in` is per kg (FullFlow energy in =
  `mass_flow_in * total_enthalpy_in`), but the header was given the enthalpy flow in W.
- `schema_version` is 2. Run archives from before this merge replay with `9b85078b`.
- Only the latest raw history is kept by default, because each is about 220 MB
  (`outputs.keep_run_history`).

Check: with the October 26 line lengths and no regen heating, the port reproduces the
committed November 3 `Outputs/run_summary.csv`: 41 of 59 metrics identical, the rest
within 0.67 %. The largest difference is COPV end pressure (1703 vs 1714 psia).

## Objective

Given a thrust history and configurable hardware, predict feed-system pressures,
temperatures, flow, inventories, and engine performance over time.

### November 9: explicit torch input topology

The requested fuel, regulated GN2 purge, and oxidizer legs now live under
`torch.network`. Fuel is FTV Inlet -> FIMV -> FICK -> FIORF -> fuel injector;
oxidizer is OTV Inlet -> OIMV -> OICK -> OIORF -> OTHX -> oxidizer injector.
Both terminate at one declared Torch Chamber, followed by one Torch Throat.
Purge is PR-REG Outlet -> FIPV -> FIPCK -> FIPORF, joining `FIORF Outlet`.
The previous TFPV/CV-3 branch and HX-1 shortcut were removed from the current
config; schema-2 archives using legacy `torch.branches` remain loadable.

OTHX is a tube specification plus a time/dT curve, with time zero at the torch
opening command as confirmed by the user. Interpolation is linear with held
endpoints; no heating is applied before opening. Missing hardware data remain
null, while inherited timing/roughness defaults are labeled placeholders.
The new topology is input-only at this stage. Its dedicated transient builder,
chamber/throat equations, and nitrogen/liquid purge junction are pending; enabling
it is rejected explicitly rather than falling back to the old main-chamber
branch approximation. `Backend/Hopper_Torch_Config.py` validates the topology
and defines the OTHX temperature interpolation helper. The hardware inventory
and parameter workbook descriptions identify the configured-but-unsolved parts.
All 37 tests pass.

### November 9: independent torch purge transient

`isolated_sims/Hopper_Torch_Purge.py` now runs the nitrogen purge independently of
the vehicle simulation. `Backend/Hopper_Torch_Purge.py` solves each restriction
in TFPV -> FIPCK -> FIPORF -> shared fuel injector for a common compressible gas
flow, then integrates chamber mass and internal energy with a bidirectional
geometric throat to the main chamber. The tee remains downstream of FIORF.
FIPV was renamed TFPV at the user's request. FIORF/FIPORF/OIORF temporarily use
8 mm2 each; the throat has no added discharge coefficient.

The selected scenario prescribes 150 psia in the torch before 2 s. At 2 s,
main backpressure steps from 14.8 to 300 psia, FIMV/OIMV receive close commands,
and TFPV receives an open command. Per-valve `command_events_s` use absolute
simulation seconds; actual openings follow full-stroke open/close travel times.
TFPV closes at 7 s as an editable example, and the run ends at 7.5 s. Each run
archives original and replayable configs, resolved settings, time histories,
station pressures, pressure ladder snapshots and transition plots.

This is explicitly a **cold-gas purge model**. Initial chamber and supply/main
reservoir temperatures are set to an assumed 293 K with constant gamma 1.4.
The regulated source is held constant (764.8 psia in this case). No combustion,
liquid displacement, closing-valve residual propellant, hot-gas chemistry,
check-valve cracking pressure, line gas storage, regulator droop or COPV
depletion is represented. Feed-valve positions are displayed without predicting
their propellant flows. OTHX heating is unused in this independent calculation.
Full torch coupling into `Hopper.py` remains pending; keep `torch.enabled=false`.

Validation compares choked/subcritical orifice flow with FullFlow, checks common
flow through each series restriction, closed/reverse supply paths, valve travel
and reversal, the exact backpressure event, finite-volume backflow and mass
conservation, and timestep refinement. With these inputs, the case completes in
seconds, gives approximately 1.286 g/s purge and 300.286 psia peak computed torch
pressure, and predicts about 15.9 mg initial main-to-torch backflow. These are
conditional cold-gas results, not a combustion or hot-purge validation.

### November 9: consolidated torch sequence, connectivity diagram and pressure margin

All independent-run commands now live in `torch.standalone_purge.sequence`: one
ordered event list changes valve targets and main pressure, releases the imposed
torch pressure, and sets the run end. Component definitions keep physical travel
times only. Old per-valve run archives still replay; mixed formats fail explicitly.
`--check-inputs` prints the resolved timeline and `torch_sequence.csv` saves it.

`Backend/Hopper_Torch_Diagram.py` generates SVG/PNG drawings and a component-edge
CSV from the loaded network. It draws the purge tee after FIORF explicitly, shows
the shared injector and chamber, and marks configured-only feeds with dashed lines.
`--diagram-only` produces an input snapshot and diagram without a transient solve.

History now includes torch-minus-main pressure, with a dedicated pressure-margin
plot. A separate fully-open steady equilibrium solves inlet flow = throat outflow,
with energy balance setting cold-gas chamber temperature to supply temperature.
The current case gives 0.28637 psi positive steady margin and an unchoked throat;
the pressure step still causes initial backflow. No minimum required pressure
margin has been supplied, so the diagnostic does not assign a design pass/fail.
Tests cover centralized events and invalid/ambiguous timelines, archive-compatible
schedules, steady flow balance, the effect of reduced throat area, and tee edges.

Purge plots now omit the pre-purge firing phase from the station ladder, sample the
valve opening/closing and full purge lifetime automatically, and include an expanded
torch-to-main ladder and PNG station histories. Both 30-second, 300/160 psia comparison
cases have been re-exported to PNGs from their saved histories without repeating solves.

### November 9 publication snapshot and interpretation

The current base config preserves the user's later fuel-injector CdA of **0.8051
mm2** and FIPORF CdA of **50 mm2**. These differ from the 0.1051 mm2 injector used
in the archived 30-second comparison. The later-input run in
`Outputs/Torch_Purge/20261110T005911_735990Z/` (UTC name, November 9 EST) reports
9.5020 g/s and 315.6931 psia torch pressure against 300 psia main, with 47.9363 g
source gas consumed in the example 2-to-7-second command window including travel.
The base timeline remains that five-second example; thirty-second cases live in
separate snapshots under `Outputs/Torch_Purge/comparison_30s_20261109T183407_098289Z/`.
Those smaller-injector cases give 1.2860 g/s and about 38.6759 g source consumption
at either 300 or 160 psia main, with steady pressure margins 0.2864 and 0.5371 psi.

A positive steady pressure margin does not by itself prove acceptable hot-gas
ingestion or thermal exposure. The discussed +/-5% peak fluctuation assumption
and +/-1%/10% sensitivity cases are exploratory; no waveform or acceptance
criterion is implemented. Frequency, transient reverse-flow mass, local mixing
and hot-gas properties remain unresolved. Source-pressure scaling and a possible
downstream purge connection were discussion only, not implemented hardware changes.
The handoff records the delivered modules, exact scenario differences, replay
instructions, archived results and model boundaries. Automated physical regression
inputs are fixed independently of the user's editable injector CdA (54 tests).

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
- [x] JSON configuration for vehicle, ordered hardware paths, engine, sequence, and solver/output settings.
- [x] Imperial tube dimensions and throat area inputs; diagram also reports throat diameter.
- [x] Compact colored schematic with component symbols and orthogonal connections.
- [x] Remove Hdot and quality panels from standard plots and their plot CSVs.
- [x] Per-run input archives, copied efficiency/map assets, and resolved derived settings.

## JSON configuration and reproducibility

`Hopper_Config.py` validates the JSON and converts explicit unit-suffixed fields
to the solver's SI inputs. `Hopper_Config_Fields.py` contains only field bindings
and unit conversions. Unknown keys and duplicate JSON keys are rejected. Optional
component defaults are materialized in the saved runnable config. Existing
hardware settings were migrated, including the user's two bottles, 3/8-inch
liquid tubes, and 100-to-600 lbf ramp. Gas-line OD/wall remain unknown; their
existing IDs are preserved instead of inventing wall geometry.

Ullage collapse remains the existing convection model: tank surface radius,
ullage heat-transfer coefficient, and liquid/ullage temperatures are explicit
inputs. There is no separate empirical collapse multiplier in this model.

Efficiency functions are declarative bounded linear functions or Pc/MR
polynomials, or CSV grids. JSON does not execute arbitrary Python expressions.
Custom Python components remain an extension point in the solver, but are not
serializable JSON hardware types.

Every driver solve creates a timestamped directory before solving and writes
`input_config.json` (original input), `run_config.json` (resolved defaults and
local asset paths), and `resolved_inputs.json` (config, derived timestep/profile/
geometry, asset hashes). The actual chamber map and efficiency CSVs are copied
into that directory and used for the run. The run config is embedded in HDF5.
Original results stay intact. Re-export with the archived config and
`--results-only`; input archives are not rewritten during re-export.

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

The JSON migration passes 18 unit tests, including unit conversions, invalid
fields/geometry/timing, relative CSV paths, and independent replayable archives.
A JSON-driven 1.5 s startup/burn/shutdown run completed and generated all 19
plots with CSVs, inventories and summary, plus config/assets in its own archive.
The new schematic was visually checked. The previous physics-specific validation
below remains applicable; this migration did not change the physical equations.

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

## Startup stall investigation

The current input review caught an initial fuel-tank pressure of 5,004 psia
against a 500 psia setpoint; the user confirmed this was a typo and it was
corrected to 500 psia. Regen outlet temperatures had also been entered as dT;
the user corrected the latest curve to +51/+44 K (344/337 K outlets at 293 K inlet).

A separate nonlinear stall was reproduced at the attempted 0.015-to-0.020 s
step, including with regen disabled. FullFlow's default SciPy iteration budget
lets a bad step run for thousands of iterations before its existing retry.
`Hopper_Solver.py` bounds each least-squares attempt and retains the same
residual checks and timestep-retry mechanism. With an 80-evaluation cap the
short regression passed the difficult step by retrying at 2.5 ms and reached
0.1 s. No timestep coarsening, relaxed tolerance, hardware changes, or installed
package edits were used. Property-table and network-build progress messages
distinguish setup cost from a nonlinear stall; total run timing is saved.

The final default is 20 evaluations per least-squares attempt, exposed as
`simulation.nonlinear_max_evaluations`. Full cold-start driver runs with the
corrected current inputs (11 s including the 10 s ramp, all plots) completed
in 439.6 s with a cap of 80 and 298.8 s with a cap of 20. The latter retained
the 5 ms step cap and 1e-6 residual tolerance, using 2,264 accepted steps and
123 retries versus 2,230 steps and 31 retries. Maximum accepted residual stayed
below 1e-6; total pressurant mass drift stayed below 9 mg in both runs.
Maximum demonstrated thrust changed from 2055.91 to 2062.31 N (0.31%). Final
LOX/fuel inventories differed by 5.42/3.03 g. Instantaneous tank-pressure
differences reached 20.7 psi because bang-bang pulse timing changed; this is
not a timestep-convergence study or proof of identical pressure histories.
The runtime target was met narrowly on this machine; other runs may vary.
All 21 unit tests pass. Benchmark artifacts are in ignored
`validation_outputs/runtime_fixed/` and `validation_outputs/runtime_fixed_20/`.

## Model boundaries

The chamber remains quasi-steady and zero-volume. Liquid feeds remain single phase.
Measured geometry, valve timing, regen curves, efficiencies, and COPV temperature
limits must be supplied by the user; example/default values are not hardware data.

# Hopper FullFlow

Transient propulsion model built with FullFlow and ThermoProp. A thrust history
sets the engine's required chamber pressure and propellant flows. The transient
model solves bottles, tank ullages, feed restrictions, valves, and engine
performance together, including pressure, temperature, flow, and inventories.

## Run and configure

```powershell
pip install -r requirements.txt
python Hopper.py
python Hopper.py --config path/to/my_config.json
```

From Git Bash, use `bash run_hopper.sh` (or `bash run_hopper.sh --config
path/to/my_config.json`). This selects a working project environment without
depending on terminal activation or an editor-generated Windows drive path.
It checks `.venv`, `.venv-validation`, then `venv` for installed dependencies.
For the existing validation environment, the direct Git Bash command is
`./.venv-validation/Scripts/python.exe Hopper.py`.

Edit **[hopper_config.json](hopper_config.json)**. `HOPPER_DATA.py` only loads and
converts it into SI values for existing solver modules; do not put inputs there.
`HOPPER_CONFIG` can also select a config for all CLI tools. Relative file paths
resolve against the selected config's directory, not the working directory.

The config groups are:

- `simulation`: solve/rebuild switches, maximum step, save interval, residual
  tolerance, Jacobian method, console output.
- `files`, `outputs`: chamber map/history paths, results root, plot selections,
  plot units, CSV export, station times, and diagram generation.
- `fluids`: propellant names, pressurant, and ambient pressure.
- `sequence`: valve demand start/ramp, step smoothing, and tail-off duration.
- `numerics`: flow floor, nozzle regularization, regen property-table resolution,
  and controller/valve integration resolution.
- `engine`, `chamber_map`: mixture ratio, throat/exit areas, efficiencies, nominal
  pressure, and equilibrium map bounds/resolution.
- `pressurization`: independent controller period and common-header volume.
- `maximum_thrust`: required sustained duration and target tracking tolerance.
- `thrust_profile`: selected ramp or named profile and interpolation mode.
- `bottles`: each bottle's volume, initial pressure/temperature, minimum-temperature
  limit, and outlet CdA. Add/remove entries to change the bottle count.
- `tanks`: initial liquid mass, volume, temperatures, initial pressure, setpoint,
  asymmetric deadbands, surface radius, ullage HTC, and initial controller state.
- `flow_paths`: ordered components in each liquid feed and pressurant leg.

Units are in field names: pressures use psia, pressure differences psi, lengths
inches, areas in² for nozzle throats/exits and mm² for CdA, temperatures K,
volumes L, masses kg, and time s. JSON accepts `true`, `false`, and `null`; comments
are not valid JSON. Unknown fields and duplicate keys are rejected.

The shipped settings preserve the current vehicle inputs: two 9 L bottles,
3/8-inch liquid tubes with 0.035-inch walls, and a 100-to-600 lbf, 10 s ramp.
Placeholder geometry and timing still require measured hardware data.

## Hardware and thermal inputs

Each flow-path entry has a `type`, unique `name`, properties, and optional
`outlet_node`. Supported types are `Line`, `GasLine`, `Orifice`, `Injector`,
`Regen`, `Valve`, `BangBangValve`, `CheckValve`, `ThrottleValve`, and
`CavitatingVenturi`. Bottle/tank types are `COPVSource` and `Tank`.
Adding an entry updates the network, station plots, and diagram automatically.
Custom Python builders need a code extension; JSON does not execute Python.

```json
{
  "type": "Line",
  "name": "Fuel inlet tube",
  "outer_diameter_in": 0.375,
  "wall_thickness_in": 0.035,
  "inner_diameter_in": null,
  "length_in": 12.0,
  "roughness_in": 0.0005905511811023622,
  "height_change_in": 0.0,
  "outlet_node": "FTV Inlet"
}
```

ID is calculated as OD minus twice wall. Alternatively specify ID and leave OD
and wall `null`; the diagram marks those missing measurements. If all three are
specified they must agree. Gas lines currently have only measured/configured IDs.
Every line shows OD/wall status, ID, and length in inches on the diagram.
`engine.throat_area_in2` is the nozzle throat input; diameter is derived and
shown in inches. Valve maximum CdAs are independent inputs, not automatically
resized when a tube changes.

Valve `open_time_s` and `close_time_s` are full-stroke travel times. Partial strokes
take proportionate time; zero is ideal instantaneous travel. A scheduled `Valve`
uses `schedule: [[times in s], [openings from 0 to 1]]` on the simulation clock.
Tank pressure control opens below setpoint minus low deadband and closes above
setpoint plus high deadband. The controller samples independently of valve travel.

Ullage collapse uses convection with `ullage_htc_w_m2_k` and surface area
`pi * radius²`, using the tank's liquid and ullage temperatures. There is no
separate empirical collapse-factor multiplier. Initial propellant is specified
by `liquid_mass_kg`; ullage volume follows from tank volume and liquid density.

A `Regen` entry has a hydraulic CdA, optional inertial passage length, and
`pc_psia_delta_t_k` pairs such as `[[0, 0], [100, 10], [300, 45]]`.
Temperature rise is linearly interpolated with held endpoints. It updates
downstream fuel temperature and liquid properties. The default rise remains zero.
The regen plot retains physical heat pickup; the removed Hdot panels were
enthalpy transport, a different quantity.

COPV temperature limits only flag a violation; delivery continues. `null`
disables a bottle's check. Bottle/header connection data and header volume are
configurable. The current 0.05 L header volume is a placeholder.

## Thrust and efficiency

Set `thrust_profile.selected` to `ramp` or a key under `profiles`.
The ramp uses `start_lbf`, `end_lbf`, `duration_s`, and `points`. Named profiles
use `time_s_thrust_n` pairs; the included `1dof` profile is preserved. Times must
be unique, start at zero, and thrusts must be positive.
`interpolation` is `linear` or `step`; physical valve travel is modeled separately.

The maximum-thrust report is the highest thrust maintained over the configured
sustained duration, within valid inventory and map bounds. If no sustained
saturation/target shortfall was demonstrated, this is a lower bound, not proof
of a global hardware maximum. The ramp must actually reach a limiting condition
to provide that evidence.

Each engine efficiency supports:

- `mode: "linear"`: `low`, `high`, `pc_low_psia`, `pc_span_psi`. Bounded linear
  interpolation with endpoint holding. A non-null `csv` overrides this function.
- `mode: "polynomial"`: `coefficients[i][j]` multiplies `Pc_psia^i * MR^j`.
  For example `[[0.8]]` is constant 0.8; `[[0.8, 0.01], [0.0001]]` is
  `0.8 + 0.01*MR + 0.0001*Pc_psia`.
- `mode: "csv"`: `csv` points to a file with `pc_psia,efficiency` or
  `pc_psia,mr,efficiency` columns. MR tables must form a complete rectangular grid.

Efficiencies must be in `(0, 1]`; polynomial outputs are validated during use.
CSV interpolation is linear/bilinear with held edges. Files in `inputs/` are
format examples, not calibration data. Changing efficiencies does not require
rebuilding the ideal equilibrium map; changing propellants or map bounds does.

## Run archives and results

Every new driver solve creates `Hopper_Results/<UTC timestamp>/` (or the configured
results root). It prints the folder before solving. Previous runs are preserved.
That folder contains:

- `input_config.json`: original input text.
- `run_config.json`: runnable config with component defaults filled in and
  paths referring to this archive's input assets.
- `resolved_inputs.json`: config plus derived throat diameter, expansion ratio,
  timestep cap, sampled thrust demand, header initialization, and asset hashes.
- `inputs/`: the actual chamber map and selected efficiency CSVs used by the run.
- FullFlow HDF5 history, also containing the run config as metadata.
- Engine targets, actual performance, residual inventory CSV, summary JSON/text,
  selected plots and optional plot CSVs, component diagram PNG/SVG.

Run snapshots are written before the transient solve, so failed attempts retain
their inputs too. A complete result also has the summary and history.

To regenerate outputs for a saved run, use its archived config:

```powershell
python Hopper.py --config Hopper_Results/<timestamp>/run_config.json --results-only
```

To make selected plots directly:

```powershell
$env:HOPPER_CONFIG = "Hopper_Results/<timestamp>/run_config.json"
python Hopper_Plots.py list
python Hopper_Plots.py plot ox_pressure_stations fuel_pressure_stations --times 1 5 --csv
python Hopper_Plots.py plot ox_ladder engine regen copv_temperatures
```

Other plots include `nodes`, `feed_nodes`, `fuel_ladder`, `pressurant_ladder`,
`stiffness`, `exit_pressure`, `ox_stations`, `fuel_stations`, `pressurization`,
`tank_states`, `startup`, `lox_pressurant_stations`, and `fuel_pressurant_stations`.
Standard node/ladder plots no longer include quality or Hdot. Raw solver tracks
remain available. Station snapshots use component names along the flow path;
gas-group intermediate pressures are allocated approximations, not separate
compressible node solves. `station_times_s: null` selects peak delivered thrust.

## Numerical limits and validation

The timestep is capped by controller and valve resolution. Current defaults
produce 5 ms steps from a requested 50 ms maximum. These settings are resolution
heuristics, not automatic temporal error estimates; check timestep convergence.
The smaller timestep and individual bottle states increase runtime.

`simulation.nonlinear_max_evaluations` limits each SciPy nonlinear solve attempt
(default 20). Difficult steps then use FullFlow's existing smaller-timestep retry
instead of spending thousands of iterations on one attempt. The residual
tolerance is unchanged; a failed step is never accepted just because its work
budget ran out. Numerical Jacobian probes are additional evaluations, so this
is not a wall-clock timeout. `Hopper_Solver.py` adapts the pinned FullFlow 3.0.0
step solver without modifying installed packages.

On the checked 11 s simulated run (10 s ramp plus sequencing), the 20-evaluation
budget completed setup, solve and all plots in 298.8 s versus 439.6 s with 80.
This is one measured run, not a runtime guarantee for other hardware/computers.
Both kept accepted residuals below 1e-6 and total pressurant mass drift below
9 mg. Smaller retries changed bang-bang pulse timing: instantaneous tank pressures
differed by up to 21 psi between the two runs, while demonstrated maximum thrust
changed 0.31% and final tank inventories by less than 6 g. Check timestep convergence
before using precise pressure pulse timing as a design result.

Startup now reports heated-liquid property-table progress and distinguishes
network construction from transient integration. Nonzero regen curves require
expensive enthalpy integrations when building this table. The driver writes
`run_timing.json` and prints total runtime including setup and output generation.
`pc_psia_delta_t_k` is a rise, not an absolute outlet temperature: with a 293 K
inlet, a 51 K rise gives a 344 K outlet. Values below the curve's first pressure
hold that endpoint, including startup; add a zero-pressure point if that is
the desired model. Positive tube `height_change_in` means an upward run.

The chamber remains quasi-steady and zero-volume. Liquid feeds are single-phase;
ignition, wall thermal inertia, two-phase depletion, and trapped liquid in
zero-volume junctions are not modeled. Regen heating does not regenerate the
equilibrium combustion map. Negative inventory invalidates the reported run.
Residual propellant is the remaining tank inventory.

```powershell
python -m unittest discover -s tests -v
python tests/run_feature_validation.py
```

See [DESIGN.md](DESIGN.md) for implementation decisions and progress. Core files:
`Hopper_Config.py` loads inputs, `Hopper_Elements.py` builds hardware,
`Hopper_Prop_System.py` solves the network, `Hopper_Components.py` models engine
performance, `Hopper_Chamber_Map.py` builds the gas map, `Hopper_Thrust_Targeter.py`
solves demands, and `Hopper_Export.py`, `Hopper_Report.py`, `Hopper_Plots.py`,
`Hopper_Diagram.py`, and `Hopper_Run_Config.py` produce results and archives.
`comparison/` retains historical comparisons with the previous vehicle model.

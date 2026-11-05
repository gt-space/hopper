# Hopper Vehicle Fluids (FullFlow)

Transient model of the Hopper vehicle fluids and engine P&ID, built with
[FullFlow](https://github.com/saakethramoju/FullFlow). It replaces
`hopper/sims/general_fluid_network/vehicle_sim.py` (the "old model" below) for
the new P&ID. Its layout follows the FullFlow author's Vespula and Raptor
examples: a data file, a chamber-map builder, a prop-system function, an export
module and a driver.

## Running

```bash
pip install -r requirements.txt
python Hopper.py
```

**Every input is in [`HOPPER_DATA.py`](HOPPER_DATA.py)**, the one file to
edit. No other file holds a vehicle number, and values marked `GUESS` have not
been measured yet. It is laid out in sections:
- run settings: time step, tolerances, output and plot options;
- thrust curves, with a selector (`PROP_THRUST_CURVE_NAME`);
- fluids: one name per propellant (`OXIDIZER = "LOX"`, `FUEL = "RP-1"`), used
  for both the feed-line liquid and the combustion map;
- COPVs;
- tanks, each with its own bang-bang controller;
- the pressurant legs and feed lines, written as P&ID element lists (below);
- engine: MR, throat, nozzle, efficiencies and the chamber-map grid;
- start-up sequence, plot-only values and solver numerics.

### The P&ID as element lists

Each flow path is an ordered list of elements, from tank (or COPV) to chamber
(or ullage):

```python
LOX_FEED = [
    Line("LOX Tank to OTV Line", inner_diameter=LOX_TUBE_ID, length=1.0, outlet_node="OTV Inlet"),
    ThrottleValve("OTV", max_cda=..., outlet_node="OTV Outlet"),
    Line("OTV to LOX Injector Line", inner_diameter=LOX_TUBE_ID, length=0.2032),
    Injector("OINJ", cda=15.3 * milli**2, length=INJECTOR_PASSAGE_LENGTH),
]
```

`Hopper_Prop_System.py` builds each element from native FullFlow components
and joins neighbours with a FullFlow junction. To add hardware, add an element
to the list; the network, tracks and plots follow automatically. The element
types are in [`Hopper_Elements.py`](Hopper_Elements.py):

| Element | Builds (FullFlow) |
|---|---|
| `Line` | `DarcyWeisbach` with flow inertia + `Churchill` friction |
| `Orifice`, `Injector` | `DischargeCoefficient` (optional passage length = inertia) |
| `Valve` | `DischargeCoefficient` with its opening driven by a `Sequence` schedule |
| `ThrottleValve` | `DischargeCoefficient` whose CdA a `Balance` sets to hit the thrust-profile flow, capped at fully open |
| `CavitatingVenturi` | `CavitatingVenturi`, with the propellant's vapour pressure |
| `GasLine` | equivalent CdA in the leg's `SeriesCdA`; `inertia=True` builds a `FlowTube` instead |
| `BangBangValve`, `CheckValve`, gas `Orifice` | members of the leg's `SeriesCdA` feeding one `CompressibleOrifice` |
| `Custom` | any FullFlow component (or your own `Component` subclass) you build |

Liquid properties come from ThermoProp `Propellant` objects, evaluated once
per tank at its set pressure and liquid temperature: the liquids are treated
as incompressible. ThermoProp's RP-1 is 812 kg/m³, against the 752 kg/m³
n-Dodecane surrogate vehicle_sim.py used, so fuel-side pressure drops are about
8 % higher than in the old model.

The defaults fly the full 1dof thrust profile (29.7 s), using linear flow-demand
interpolation. `PROP_DT` is the requested maximum step. Valve/controller resolution
currently reduces it to 0.005 s; the selected integration step is printed.
`PROP_PROFILE_INTERPOLATION = "step"` instead holds demands between profile points.
Runtime depends strongly on the resolved step and hardware settings. The chamber map (`Hopper_Chamber_Map.h5`) is
included, and is rebuilt automatically if missing or when
`PROP_GENERATE_CHAMBER_MAP = True` (~2 min). Rebuild it after changing the
propellants or the map grid.

Outputs go to `Hopper_Results/`:
- `engine_targets.csv` and `actual_performance.csv`, in vehicle_sim.py's formats
  (plus a target-thrust column);
- the plots named in `PROP_PLOTS` (default `"all"`, including station and thermal plots), and
  with `PROP_PLOT_CSV = True` the numbers behind each one.

The full FullFlow history goes to `Hopper_Prop_System.h5` (~340 MB). Keep it to make more plots later without re-running.

## Plots

`Hopper_Plots.py` makes every graph in the old sim's output deck
(`outputs needed.pptx`) from the saved results, plus six more. It can also
plot any node, element, tracked quantity or component attribute on request.
`Hopper.py` makes the `PROP_PLOTS` set at the end of each run.

### The deck's graphs

| Deck slide | Old sim figure | Plot name |
|---|---|---|
| 2, 9 | Network Simulation: Nodes (pressure, temperature, mass, density, quality, fill level) | `nodes` |
| 3, 11 | Ox Feed System pressure ladder (mdot, dP, CdA, Hdot, Q, State) | `ox_ladder` |
| 4, 11 | Fuel Feed System pressure ladder | `fuel_ladder` |
| 5 | Pressurization System pressure ladder | `pressurant_ladder` |
| 6, 10 | Engine Transient Performance (Pc, thrust vs target, MR) | `engine` |
| 6, 10 | Injector Stiffness | `stiffness` |
| 7 | Nozzle Exit Pressure vs Throttle | `exit_pressure` |

`plot deck` makes all seven. Slides 9-11 show the same figures for a thrust
ramp that this model does not include; any curve added to `THRUST_CURVES` in
`HOPPER_DATA.py` gets the same plots.

The other named plots are:
- `feed_nodes`: the node grid for every station along both feed lines;
- `ox_stations` and `fuel_stations`: station pressures, plus the OTV or FTV
  CdA as a percentage of fully open;
- `pressurization`: COPV and tank pressures, and the bang-bang states;
- `tank_states`: gas temperatures and remaining masses;
- `startup`: the main-valve opening over the first 2 s.

### On request

```bash
python Hopper_Plots.py list                                     # every plot, node, element, track and component
python Hopper_Plots.py plot deck                                # the deck's seven graphs
python Hopper_Plots.py plot ox_ladder engine --profile-time     # named plots, timed like the old sim
python Hopper_Plots.py nodes "LOX Tank" COPV --csv              # node grid for any nodes or node groups
python Hopper_Plots.py ladder OTV OINJ "Regen Jacket" --tmin 5 --tmax 8   # ladder for any elements, zoomed
python Hopper_Plots.py track "Chamber Pressure [psia]" "Thrust [N]"       # any tracked quantity
python Hopper_Plots.py component OTV mass_flow upstream_pressure          # any stored component attribute
python Hopper_Plots.py plot exit_pressure --mr 1.25             # exit-pressure sweep at another MR
```

Options for every command:

| Option | Effect |
|---|---|
| `--file` | Results file to read (default `Hopper_Prop_System.h5`). |
| `--out` | Output folder (default `Hopper_Results`). |
| `--units E` / `--units SI` | psi and deg F, as the old sim (default), or Pa and K. |
| `--csv` | Also write each plot's numbers to `<plot>.csv`. |
| `--profile-time` | Time from the start of the thrust profile, the old sim's clock. This model's clock starts 0.5 s earlier, with the main-valve ramp. |
| `--tmin`, `--tmax` | Time window. |
| `--name` | File name for a custom plot. |
| `--show` | Also open the plot windows. |

Node groups are `tanks`, `ox_feed`, `fuel_feed` and `all`. Element
groups are `ox`, `fuel` and `pressurant`. Bad names stop with
a one-line error that lists the valid ones. From Python:

```python
from Hopper_Plots import HopperRun, make_plots
make_plots(HopperRun("Hopper_Prop_System.h5"), ["nodes", "ox_ladder"], "Hopper_Results", csv_out=True)
```

### What the plotted quantities are

The model tracks every element's flow, pressures, CdA, state and inlet
enthalpy, and every junction's pressure. The plots read those tracks with
FullPlot, and their node and element catalogue is generated from the element
lists in `HOPPER_DATA.py`, so added elements show up automatically.
Properties are ThermoProp's, as in the model.

- **Nodes.** Tees, valve outlets and manifolds have no volume, so their mass
  is 0. Their fill level is 1 when liquid-full and 0 for gas nodes. Quality
  is blank, with a note, because every state in this model is single-phase.
  The old sim drew those states as 0.
- **Lines.** A line's CdA is its friction-only equivalent, A / sqrt(f L / D),
  using the solved friction factor. This is how vehicle_sim.py's `Line`
  defines it.
- **Pressurant legs.** Each run of gas restrictions is one FullFlow
  `CompressibleOrifice` at their `SeriesCdA`. The ladder splits its dP between
  them with the old sim's `Series` rule, dP_i = dP_group (CdA_series / CdA_i)^2.
  While the solenoid is shut, it carries the whole dP.
- **Hdot** is mdot times the enthalpy of the entering fluid, on ThermoProp's
  reference. For the propellants that includes the heat of formation, so its
  values differ from the old sim's, which used CoolProp's reference.
- **State** is the old sim's definition:
  - 1 for fixed restrictions;
  - the open/closed command for the bang-bang solenoids;
  - for OTV/FTV, the commanded flow / full-scale flow, as vehicle_sim.py's
    `mdot`-mode `ThrottleValve` defined it;
  - for the tank-to-chamber entry, 1 while the feed is commanded open.

  The full-scale flows are 0.777 and 0.388 kg/s (`OTV_FULL_SCALE_MASS_FLOW`
  and `FTV_FULL_SCALE_MASS_FLOW`). The deck's own runs appear to have used
  ~0.756 kg/s for OTV, so its OTV State reads ~3 % higher for the same flow.
  How far OTV/FTV are actually open is in `ox_stations` / `fuel_stations`:
  100 % there means the valve is saturated.
- **Exit pressure** is a steady sweep of the same engine block the transient
  uses. Throttle is Pc / `NOMINAL_CHAMBER_PRESSURE` (300 psia). The sweep
  runs at `MIXTURE_RATIO` (1.8) unless `--mr` is given, and the run's own
  operating points are drawn over it. It reproduces the deck's figure to
  within 1 % (15.81 against 15.88 psia at 100 %).

## Comparison with vehicle_sim.py

`comparison/` reruns the old model and compares it with this one on the 1dof
profile.

```bash
# 1. Old model: needs the general_fluid_network requirements.txt environment
python comparison/run_vehicle_sim.py --gfn-dir <hopper>/sims/general_fluid_network --output-dir comparison/vehicle_sim_outputs
# 2. This model
python Hopper.py
# 3. Overlay plots and the difference table
python comparison/compare_models.py --old-history comparison/vehicle_sim_outputs/vehicle_sim_history.csv --new-results Hopper_Prop_System.h5 --output-dir comparison/results
```

`run_vehicle_sim.py` runs vehicle_sim.py unchanged on disk, with two in-memory
changes:
- the selected profile is swapped from `thrust_curve_throttle_test` to
  `thrust_curve_1dof`;
- matplotlib uses a non-interactive backend, so figures are saved instead of
  shown.

It also dumps every history the old model keeps to `vehicle_sim_history.csv`.

## Configuring the new capabilities

The working design and task checklist are in [DESIGN.md](DESIGN.md).

### Flight profile or ramp

Set `PROP_THRUST_CURVE_NAME` to `"1dof"` or `"ramp"`. For a ramp, set
`RAMP_START_THRUST`, `RAMP_END_THRUST` (N; multiply lbf by `LBF`),
`RAMP_DURATION` (s), and `RAMP_POINTS`. The targeter solves each supplied thrust
point; the transient interpolates the resulting flow demands.

After a new run, `run_summary.txt` and `run_summary.json` report maximum
demonstrated sustained thrust. The value is the highest **minimum** thrust
maintained for `MAX_THRUST_HOLD_TIME`, considering only valid profile samples.
This deliberately rejects brief startup spikes. The report distinguishes a run
that demonstrates sustained valve saturation plus thrust shortfall from a run
that stops below the feed-system capacity. Neither result is a universal steady
maximum independent of initial pressure, ramp rate, mixture ratio, or depletion.
The configured ramp endpoint is a demand, not a guaranteed capability.

### Valve travel and tank deadband

`Valve`, `ThrottleValve`, `BangBangValve`, and `CheckValve` accept `open_time` and
`close_time` in seconds. These are full-stroke linear travel times; a partial
stroke takes proportionately less time. Zero means instantaneous. The example
bang-bang values are 0.05 s and require replacement with measured hardware data.
There is no extra solenoid electrical delay or valve-specific nonlinear CdA
versus position curve in this model.

`BANG_BANG_SAMPLE_PERIOD` sets the pressure-controller sample clock independently.
`LOX_TANK_SET_PRESSURE` / `FUEL_TANK_SET_PRESSURE` set the control targets. Each
has separate `...DEADBAND_LOW` and `...DEADBAND_HIGH` values in Pa. A 5 psi lower
and upper offset gives a total 10 psi band. Opening/closing travel follows the
controller command rather than switching the hydraulic area instantly.

The solver caps the requested `PROP_DT` to at most 1/10 of the shortest nonzero
valve travel and 1/2 the controller sample period. Controller sample instants
are supplied as solver breakpoints. This is an initial resolution rule;
check time-step convergence for the quantities you care about. `PROP_SAVE_DT`
separately controls saved resolution: use e.g. 0.005 s when studying 50 ms valves.
`THROTTLE_STEP_TIME` smooths the demand in step mode; it is not actuator travel.
The main-valve startup/shutdown demand ramps remain independently configurable.

### Regen temperature rise

Edit `FUEL_REGEN_DT_VS_PC`, a list of `(Pc in Pa, dT in K)` pairs. For example:

```python
FUEL_REGEN_DT_VS_PC = [(0.0, 0.0), (100 * psi, 10.0), (300 * psi, 45.0)]
```

Values are linearly interpolated and held at the endpoints. The default curve
has zero rise because no measured curve has been supplied. The `Regen` element
retains its hydraulic CdA and passage inertia. It calculates outlet temperature
and updates downstream fuel properties for the injector and subsequent elements.
The `regen` plot shows inlet/outlet temperature and prescribed heat pickup.

This is a prescribed, single-phase heating model. It has no wall thermal inertia
or reverse-flow heat transport; keep inputs within liquid-property validity.
The equilibrium chamber map still uses its original reactant temperatures, so
the prescribed rise does not automatically change equilibrium gas properties.

### COPVs and residual propellant

`COPVS` contains individual `COPVSource` objects. Each has volume, initial
pressure/temperature, a bottle-to-header `outlet_cda`, and optional
`min_temperature` (K). `COPV_MIN_TEMPERATURES` supplies one limit per default bottle;
`None` disables that bottle's limit. Bottle connection CdAs are placeholders and
should be replaced by actual connection data. Individual bottles feed a common
adiabatic header with gas storage; backflow into bottles is blocked.
`COPV_HEADER_VOLUME` sets its internal volume. The 0.05 L default is a placeholder
for measured manifold/plumbing volume. A positive volume keeps pressure defined
when both pressurant valves close and stores the gas trapped in that header.

Minimum temperature is a **reported operating limit**, not an artificial
temperature clamp or shutoff. The run continues after crossing it.
`copv_temperatures` plots all bottles and their limits. Aggregate legacy COPV
pressure/temperature tracks now describe the common header; aggregate mass is
the sum of all bottle masses.

`propellant_residuals.csv` reports initial, consumed, profile-end, and final liquid
mass plus percentage remaining. Residuals refer to tank inventory; the model
does not count trapped liquid in zero-volume feed junctions. Negative inventory
invalidates the run and is reported; two-phase depletion is not simulated.

### Efficiency functions or CSV tables

Edit `CSTAR_EFFICIENCY_FUNCTION(pc, mr)` and `CF_EFFICIENCY_FUNCTION(pc, mr)`;
functions receive chamber pressure in **Pa** and O/F mass ratio. Alternatively set
`CSTAR_EFFICIENCY_CSV` / `CF_EFFICIENCY_CSV` to paths relative to this repo.
A selected CSV takes precedence over the corresponding function.

CSV formats:

```csv
pc_psia,efficiency
100,0.80
300,0.85
```

or `pc_psia,mr,efficiency` for a complete rectangular Pc/MR grid. Interpolation is
linear/bilinear, with endpoint holding. Missing grid cells, duplicate points,
nonfinite values, and efficiencies outside `(0, 1]` are rejected. Example files
in `inputs/` describe the format; their values are not measured calibration.

The default functions now interpolate LOW to HIGH and hold their endpoints.
This corrects the former `max(1, ...)` formula, which held HIGH at low Pc and
increased without a cap. Expect results to differ from old saved comparisons.
Changing efficiencies does not require rebuilding the equilibrium gas map.

### Station snapshots and component diagram

```powershell
python Hopper_Plots.py plot ox_pressure_stations fuel_pressure_stations --times 1 5 10 --csv
python Hopper_Plots.py plot lox_pressurant_stations fuel_pressurant_stations --times 5
python Hopper_Plots.py plot regen copv_temperatures
python Hopper_Diagram.py
```

Station x labels are the source followed by each component's outlet. Times are
simulation seconds unless `--profile-time` is used. `STATION_PLOT_TIMES = None`
selects the time of peak delivered thrust. Pressure drops within grouped gas
restrictions are allocated by their instantaneous effective CdAs; these are not
separate compressible node solves.

`component_diagram.png` and `.svg` show every configured source and element.
Liquid tube boxes include OD/wall when configured, and always show ID and length.
Gas lines configured only by ID show ID; OD/wall are not invented. `Hopper.py`
also generates this diagram when `PROP_COMPONENT_DIAGRAM` is true.

New runs embed profile timing, limits, initial inventories, and throttle areas in
the HDF5 file for reliable reports. Older histories can still make the existing
plots and station snapshots, but require rerunning to obtain the new summary.

### Validation

```powershell
python -m unittest discover -s tests -v
python tests/run_feature_validation.py
```

These checks cover timing/retries, sampling, ramp generation, CSV interpolation
and validation, and sustained-thrust reporting. The feature-validation command
runs a short coupled burn with temporary in-process overrides and checks
inventory conservation, regen outlet properties, valve rates, limits, report
exports, and every named plot. It writes to `validation_outputs/` and does not
overwrite the normal saved run.

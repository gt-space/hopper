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

### Folder layout

| Where | What | Edit it? |
|---|---|---|
| `HOPPER_DATA.py` | every input: components, properties, routing, sequence, Monte Carlo bands | yes, this is the one file to edit |
| `thrust_curves/` | thrust-curve CSVs | add new curves here |
| `Hopper.py`, `Hopper_MonteCarlo.py` | the two scripts you run | no |
| `ASSUMPTIONS.txt` | every modelling assumption and its source | keep it current |
| `Backend/` | the model and its tools (element classes, network builder, engine, targeter, plots, report, workbook) | only to add a new element type or output |
| `isolated_sims/` | standalone calculators: tank vent time, LOX drain time | no |
| `Old_sims_comparison/` | reruns vehicle_sim.py and compares it with this model | no |
| `Outputs/` | everything a run writes (CSVs, plots, workbook, `Hopper_Prop_System.h5`, `MonteCarlo/`) | no, regenerated |

Every path is anchored to this folder, so scripts work from any working
directory. Scripts inside the subfolders run directly, e.g.
`python Backend/Hopper_Plots.py list` or `python isolated_sims/Hopper_Vent.py`.

**Every input is in [`HOPPER_DATA.py`](HOPPER_DATA.py)**, the one file to
edit. No other file holds a vehicle number, and values marked `GUESS` have not
been measured yet. It is laid out in sections:
- run settings: time step, tolerances, output and plot options;
- thrust curves, as CSV files in `thrust_curves/`, with a selector
  (`PROP_THRUST_CURVE_NAME`);
- fluids: one name per propellant (`OXIDIZER = "LOX"`, `FUEL = "RP-1"`), used
  for both the feed-line liquid and the combustion map;
- COPVs;
- tanks, each with its own bang-bang controller and ullage collapse factor;
- the pressurant legs and feed lines, written as P&ID element lists (below);
- side paths: every vent, relief, purge and drain on the P&ID, plus the
  regulator, as `Branch` lists;
- engine: MR, throat, nozzle, efficiencies and the chamber-map grid;
- start-up sequence, Monte Carlo dispersions, report checks and solver numerics.

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
| `Branch` | a side path off a node (COPV, ullage, tank liquid, regulator outlet) to ambient, the chamber or a feed junction: one `SeriesCdA` into one `CompressibleOrifice` (gas) or `DischargeCoefficient` (liquid) |
| `ReliefValve` | branch restriction opening linearly from crack to full-open dP |
| `Regulator` | ideal regulator off a gas node; its outlet (`"PR-REG Outlet"`) feeds branches |
| `Valve(..., schedule=CLOSED)` | a valve shut for the whole run (vents, purges, fill/drain) |

Side paths look like this:

```python
BRANCHES = [
    Branch("OVNT Vent", "LOX Ullage", [Valve("OVNT", cda=9.4 * milli**2, schedule=CLOSED)]),
    Branch("ORV Relief", "LOX Ullage", [ReliefValve("ORV", cda=30.8 * milli**2, crack_pressure=550 * psi)]),
    Branch("LOX Purge", "PR-REG Outlet", [Valve("OPV", ...), CheckValve("CV-5", ...)], sink="LOX Injector Manifold"),
]
```

Give a vent or purge a `(times, openings)` schedule to open it during a run.
A branch can also start at a feed junction and, with `into_sink=True`,
deliver its liquid into the main chamber: the torch igniter is built this
way (OIMV taps above OTV through HX-1, FIMV above FTV, firing on
`TORCH_SCHEDULE`), and its propellant counts in the chamber flow and MR.
P&ID items with no flow path in flight (quick disconnects, diffuser, the
torch spark, sensors) are listed in `NOT_MODELLED` so they still appear in
the component table and schematic. Every assumption is in `ASSUMPTIONS.txt`.

**Ullage collapse** is a collapse factor per tank (FullFlow's Vespula
pattern): the COPV supplies the full pressurant flow and the ullage keeps
`1 / collapse_factor` of it. LOX uses 2.0 and fuel 1.1 (assumed values;
Luca is to supply conservative and Vespula-heritage numbers).

**Thrust curves** are CSV files (`thrust_curves/1dof.csv`) with a header row.
`Hopper_Inputs.load_thrust_curve` finds the time and thrust columns, reads
units from the header (s or ms; N, kN or lbf) and resamples dense
trajectory exports onto `THRUST_CURVE_MIN_SPACING`. To fly a new export, drop
the file in `thrust_curves/`, add it to `THRUST_CURVES` and select it.

Liquid properties come from ThermoProp `Propellant` objects, evaluated once
per tank at its set pressure and liquid temperature: the liquids are treated
as incompressible. ThermoProp's RP-1 is 812 kg/m³, against the 752 kg/m³
n-Dodecane surrogate vehicle_sim.py used, so fuel-side pressure drops are about
8 % higher than in the old model.

The defaults fly the full 1dof thrust profile (29.7 s) at a 0.05 s time step,
with the throttle valves stepping at each profile point the way vehicle_sim.py
does (`PROP_PROFILE_INTERPOLATION = "step"`; `"linear"` interpolates instead).
A run takes about 1 min. The chamber map (`Hopper_Chamber_Map.h5`) is
included, and is rebuilt automatically if missing or when
`PROP_GENERATE_CHAMBER_MAP = True` (~2 min). Rebuild it after changing the
propellants or the map grid.

Outputs go to `Outputs/`:
- `engine_targets.csv` and `actual_performance.csv`, in vehicle_sim.py's formats
  (plus a target-thrust column);
- the plots named in `PROP_PLOTS` (default `"all"`, the 13 plots below), and
  with `PROP_PLOT_CSV = True` the numbers behind each one;
- with `PROP_REPORT = True` (default), the run report:
  - `run_summary.csv`: burn time, total impulse, propellant residuals (kg and
    % left at shutdown), COPV pressure and gas used, tank pressure min/max,
    injector stiffness, relief venting, thrust error at each profile point,
    and whether every thrust level is reachable (OTV/FTV below 100 % of max CdA);
  - `cda_table.csv`: per thrust-profile point, target and achieved thrust, Pc,
    MR, flows, OTV/FTV CdA and % of max, injector dP, tank and COPV pressure;
  - `vent_times.csv`: each tank vent's blow-down time, full and empty tank,
    against `VENT_TIME_REQUIREMENTS` (`Hopper_Vent.py`);
  - `lox_drain_times.csv`: time to drain the loaded LOX through OFD (the same
    Jaksa valve as OVNT, plus any `LOX_DRAIN_EXTRA_ELEMENTS` such as OQD) with
    the ullage held at each pressure in `LOX_DRAIN_ULLAGE_PRESSURES`, from a
    small FullFlow transient with the liquid head (`Hopper_Drain.py`);
  - `component_table.csv` / `.md`: every P&ID item with type, CdA, set points
    and source (`Hopper_Hardware.py`);
  - `schematic.png`: the simulated network, drawn from the element lists
    (`Hopper_Schematic.py`);
  - `Hopper_Parameters.xlsx`: every input with status (GUESS / EDIT), owner
    team, every component, the modelling assumptions and the open items
    (`Hopper_Parameters.py`).

Each of these also runs on its own, e.g. `python Backend/Hopper_Report.py` on the
last results file, or `python Backend/Hopper_Schematic.py`.

The full FullFlow history goes to `Hopper_Prop_System.h5` (~340 MB). Keep it to make more plots later without re-running.

## Monte Carlo

```bash
python Hopper_MonteCarlo.py        # MC_RUNS runs, MC_WORKERS at a time
```

Each run is the full transient with inputs drawn from `MC_DISPERSIONS` in
`HOPPER_DATA.py`: groups of inputs, each with a percent (or absolute) error
band, drawn uniformly. `MC_ALL_CDAS` stands for every CdA that carries flow,
so an element added to a list is varied automatically. A group is "each"
(every input drawn separately, e.g. the CdAs) or "shared" (one draw moves the
whole group, e.g. the c* efficiency curve). Current bands: CdA +-10 %,
initial tank and COPV pressure +-2 %, loads +-1 %, collapse factor +-0.3,
c* +-2 %, Cf +-1 % (all listed in `ASSUMPTIONS.txt`, section 14).
The seed makes the draws repeatable and run 0 is nominal. With
`MC_THRUST_CONTROL = "thrust"` the thrust targets are re-solved with each
run's engine (the controller holds thrust); `"flow"` flies nominal flows.
Results in `Outputs/MonteCarlo/`: every run's inputs and summary metrics
(`monte_carlo_runs.csv`), the statistics of each metric, i.e. the error bars
(`monte_carlo_stats.csv`: mean, std, 5/50/95 %), and time-history envelopes
(`monte_carlo_envelopes.png`).

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
python Backend/Hopper_Plots.py list                                     # every plot, node, element, track and component
python Backend/Hopper_Plots.py plot deck                                # the deck's seven graphs
python Backend/Hopper_Plots.py plot ox_ladder engine --profile-time     # named plots, timed like the old sim
python Backend/Hopper_Plots.py nodes "LOX Tank" COPV --csv              # node grid for any nodes or node groups
python Backend/Hopper_Plots.py ladder OTV OINJ "Regen Jacket" --tmin 5 --tmax 8   # ladder for any elements, zoomed
python Backend/Hopper_Plots.py track "Chamber Pressure [psia]" "Thrust [N]"       # any tracked quantity
python Backend/Hopper_Plots.py component OTV mass_flow upstream_pressure          # any stored component attribute
python Backend/Hopper_Plots.py plot exit_pressure --mr 1.25             # exit-pressure sweep at another MR
```

Options for every command:

| Option | Effect |
|---|---|
| `--file` | Results file to read (default `Hopper_Prop_System.h5`). |
| `--out` | Output folder (default `Outputs`). |
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
from Backend.Hopper_Plots import HopperRun, make_plots
make_plots(HopperRun("Outputs/Hopper_Prop_System.h5"), ["nodes", "ox_ladder"], "Outputs", csv_out=True)
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

`Old_sims_comparison/` reruns the old model and compares it with this one on the 1dof
profile.

```bash
# 1. Old model: needs the general_fluid_network requirements.txt environment
python Old_sims_comparison/run_vehicle_sim.py --gfn-dir <hopper>/sims/general_fluid_network --output-dir Outputs/Old_sims_comparison/vehicle_sim_outputs
# 2. This model
python Hopper.py
# 3. Overlay plots and the difference table
python Old_sims_comparison/compare_models.py --old-history Outputs/Old_sims_comparison/vehicle_sim_outputs/vehicle_sim_history.csv --new-results Outputs/Hopper_Prop_System.h5 --output-dir Outputs/Old_sims_comparison/results
```

`run_vehicle_sim.py` runs vehicle_sim.py unchanged on disk, with two in-memory
changes:
- the selected profile is swapped from `thrust_curve_throttle_test` to
  `thrust_curve_1dof`;
- matplotlib uses a non-interactive backend, so figures are saved instead of
  shown.

It also dumps every history the old model keeps to `vehicle_sim_history.csv`.

Transient model of the Hopper vehicle fluids and engine P&ID, built with FullFlow.

## Running

```bash
pip install -r requirements.txt
python Hopper.py
python Hopper.py --config path/to/my_config.json   # another case
python Hopper.py --results-only                    # redo every output from Outputs/Hopper_Prop_System.h5
```

From Git Bash, `bash run_hopper.sh` finds a working project environment (`.venv`,
`.venv-validation`, `venv`) and forwards the same arguments.

### Folder layout

| Where | What | Edit it? |
|---|---|---|
| `hopper_config.json` | every input: components, properties, routing, sequence, Monte Carlo bands | yes, this is the one file to edit |
| `thrust_curves/` | thrust-curve CSVs | add curves here |
| `Hopper.py`, `Hopper_MonteCarlo.py` | the two scripts you run | no |
| `HOPPER_DATA.py` | loads the JSON into the Python names the model uses, in SI units | only when adding a new input field |
| `ASSUMPTIONS.txt` | key inputs and modelling assumptions | keep it current |
| `Backend/` | the model and its tools (config loader, element classes, network builder, engine, targeter, plots, report, workbook) | only to add a new element type or output |
| `isolated_sims/` | standalone calculators: tank vent time, LOX drain time, tank pressure sweep, purge flow check | no |
| `Old_sims_comparison/` | reruns vehicle_sim.py and compares it with this model | no |
| `Outputs/` | the latest run's results for viewing (CSVs, plots, workbook, `Hopper_Prop_System.h5`, `MonteCarlo/`) | no |
| `Hopper_Results/` | one folder per solve with its exact inputs and results | no |
| `tests/` | unit tests and a short coupled feature check | no |
| `DESIGN.md`, `HANDOFF.md` | design record and the 2026-11-05 JSON migration handoff | no |

### Inputs

Unit suffixes in the field names are authoritative: psia absolute, psi differences, psig
gauge, inches, in² for the throat and exit, mm² for CdA, K, L, kg, s. Unknown fields and
duplicate keys are rejected. Relative paths are relative to the config file.

JSON has no comments, so tags go in a `"note"` on a component, or in the `notes` section by
JSON path (`"engine.mixture_ratio": "EDIT target O/F ..."`). `EDIT` = expected to change,
`GUESS` = nobody has confirmed it. Notes feed the component table and the workbook.

The sections are:

- `simulation`, `files`, `outputs`, `report`: step, tolerances, paths, plot and report options;
- `fluids`: one name per propellant, used for both the feed-line liquid and the combustion map;
- `sequence`: main-valve start and ramp, thrust slew limit, tail-off;
- `engine`, `chamber_map`: MR, throat, nozzle, efficiencies and the chamber-map grid;
- `pressurization`: bang-bang sample period and the COPV model (below);
- `thrust_profile`: curves as CSV files in `thrust_curves/` or inline, plus a parameterized
  ramp, with a selector (`selected`);
- `bottles`: each COPV;
- `tanks`: each with its own bang-bang band and ullage collapse factor;
- `flow_paths`: the pressurant legs and feed lines, as P&ID element lists (below);
- `regulators`, `side_paths`, `torch`: every vent, relief, purge and drain on the P&ID;
- `vent`, `lox_drain`, `tank_sweep`: inputs of the isolated calculators;
- `not_modelled`: P&ID items with no flow path, for the component table and schematic;
- `monte_carlo`: runs, seed, dispersions.

### The P&ID as element lists

Each flow path is an ordered list of elements, for example:

```json
"oxidizer_feed": [
  {"type": "Line", "name": "LOX Tank to OTV Line", "outer_diameter_in": 0.375, "wall_thickness_in": 0.01,
   "length_in": 48, "outlet_node": "OTV Inlet", "note": "EDIT 3/8 in x 0.01 in wall"},
  {"type": "ThrottleValve", "name": "OTV", "max_cda_mm2": 50.0, "full_scale_mass_flow_kg_s": 0.777, "outlet_node": "OTV Outlet"},
  {"type": "Line", "name": "OTV to LOX Injector Line", "outer_diameter_in": 0.375, "wall_thickness_in": 0.01, "length_in": 10},
  {"type": "Injector", "name": "OINJ", "cda_mm2": 16.6, "length_in": 0.3937}
]
```

`Backend/Hopper_Prop_System.py` builds each element from native FullFlow components and
joins neighbours with a FullFlow junction. To add hardware, add an element to the list; the
network, tracks and plots follow automatically. The element types are in
[`Backend/Hopper_Elements.py`](Backend/Hopper_Elements.py):

| `type` | Built from |
|---|---|
| `Line` | `DarcyWeisbach` with flow inertia + `Churchill` friction. Give OD and wall, or ID alone |
| `Orifice`, `Injector` | `DischargeCoefficient` (optional passage length = inertia) |
| `Regen` | an `Orifice` plus a fuel temperature rise vs Pc (`pc_psia_delta_t_k`); everything downstream uses the heated properties |
| `Valve` | `DischargeCoefficient` with its opening driven by a schedule: `[[times], [openings]]`, `"closed"` (shut all run) or `"torch"` |
| `ThrottleValve` | `DischargeCoefficient` whose CdA a `Balance` sets to hit the thrust-profile flow, capped at fully open |
| `CavitatingVenturi` | `CavitatingVenturi`, with the propellant's vapour pressure |
| `GasLine` | equivalent CdA in the leg's `SeriesCdA`; `"inertia": true` builds a `FlowTube` instead |
| `BangBangValve`, `CheckValve`, gas `Orifice` | members of the leg's `SeriesCdA` feeding one `CompressibleOrifice` |
| `Branch` | a side path off a node (COPV, ullage, tank liquid, regulator outlet) to ambient, the chamber or a feed junction: one `SeriesCdA` into one `CompressibleOrifice` (gas) or `DischargeCoefficient` (liquid) |
| `ReliefValve` | branch restriction opening linearly from crack to full-open dP |
| `Regulator` | ideal regulator off a gas node; its outlet (`"PR-REG Outlet"`) feeds branches |

Valves and valve-like elements also take `open_time_s` / `close_time_s`, full-stroke travel
times (0 = instant; partial strokes take proportionate time). Nonzero travel caps the time
step at 1/10 of the shortest stroke, so it costs run time. The shipped config uses 0.

A side path looks like this:

```json
{"type": "Branch", "name": "OVNT Vent", "source": "LOX Ullage",
 "elements": [{"type": "Valve", "name": "OVNT", "cda_mm2": 9.4, "schedule": "closed"}]}
```

Give a vent or purge a `[[times], [openings]]` schedule to open it during a run. A branch
can also start at a feed junction and, with `"into_sink": true`, deliver its liquid into the
main chamber: the torch igniter is built this way (OIMV taps above OTV through HX-1, FIMV
above FTV, firing on the `torch` section's window, `torch.enabled` switches it on), and its
propellant counts in the chamber flow and MR.

### Two COPV models

`pressurization.copv_model`:

- `"lumped"` (default): all bottles as one volume. It needs bottles at the same pressure and
  temperature, and runs at the full 0.05 s step.
- `"header"`: each bottle is its own volume feeding a small common header through
  `outlet_cda_mm2`. Use it for bottles that differ. The header is stiff, so set
  `maximum_step_s` near 0.005; a run takes roughly ten times longer.

Each bottle has an optional `minimum_temperature_k`. Crossing it flags the run and does not
stop the flow.

### Efficiencies

`engine.cstar_efficiency` and `cf_efficiency` each take:

- `"mode": "linear"`: `low`, `high`, `pc_low_psia`, `pc_span_psi`, held flat outside the ramp;
- `"mode": "polynomial"`: `coefficients[i][j]` multiplies `Pc_psia^i * MR^j`;
- `"mode": "csv"`: a file with `pc_psia,efficiency` or `pc_psia,mr,efficiency` columns
  (`inputs/` has format examples, not calibration data).

Changing an efficiency does not need a new chamber map. Changing propellants, reactant
temperatures or the map grid does (`simulation.rebuild_chamber_map`).

## Outputs

Every solve creates `Hopper_Results/<UTC time>/` before it starts, with `input_config.json`,
`run_config.json` (runnable, paths into the archive), `resolved_inputs.json` (derived values
and file hashes) and `inputs/` (the chamber map, thrust curves and efficiency CSVs it used).
The results are written there, then copied to `Outputs/`, and the raw history moves to
`Outputs/Hopper_Prop_System.h5` so one large file is kept rather than one per run
(`outputs.keep_run_history` leaves a copy in every archive).

- `engine_targets.csv` and `actual_performance.csv`, in vehicle_sim.py's formats (plus a
  target-thrust column);
- `component_diagram.png` / `.svg`: compact card view of the flow paths, tube sizes in inches;
- the plots named in `outputs.plots` (default `"all"`), and with `plot_csv` the numbers
  behind each one;
- `run_summary.txt` / `.json`, `propellant_residuals.csv`: sustained demonstrated thrust,
  residuals, COPV temperature limits, validity flags;
- with `outputs.report` (default), the run report:
  - `run_summary.csv`: burn time, total impulse, propellant residuals (kg and % left at
    shutdown), COPV pressure and gas used, tank pressure min/max, injector stiffness, relief
    venting, thrust error at each profile point, and whether every thrust level is reachable
    (OTV/FTV below 100 % of max CdA);
  - `cda_table.csv`: per thrust-profile point, target and achieved thrust, Pc, MR, flows,
    OTV/FTV CdA and % of max, injector dP, tank and COPV pressure;
  - `vent_times.csv`: each tank vent's blow-down time, full and empty tank, against
    `vent.time_requirements_s` (`isolated_sims/Hopper_Vent.py`);
  - `lox_drain_times.csv`: time to drain the loaded LOX through OFD (plus any
    `lox_drain.extra_elements` such as OQD) with the ullage held at each test pressure, from a
    small FullFlow transient with the liquid head (`isolated_sims/Hopper_Drain.py`);
  - `tank_pressure_sweep.csv` / `.png`: the throttle-valve CdA each tank set pressure needs to
    fly every point of the thrust profile, against the valve and a design margin; gives the
    lowest usable tank pressure (`isolated_sims/Hopper_Tank_Pressure_Sweep.py`, steady state);
  - `purge_flow.csv` / `.png`: whether the LOX and fuel GN2 purges can match each propellant's
    peak volumetric flow, with the COPV at its end-of-run state
    (`isolated_sims/Hopper_Purge.py`);
  - `component_table.csv` / `.md`: every P&ID item with type, CdA, set points and source;
  - `schematic.png`: the simulated network, drawn from the element lists;
  - `Hopper_Parameters.xlsx`: every input with status (GUESS / EDIT), owner team, every
    component, the modelling assumptions and the open items.

Each isolated calculator also runs on its own, for example `python isolated_sims/Hopper_Drain.py`.

## Monte Carlo: WORK IN PROGRESS

```bash
python Hopper_MonteCarlo.py        # monte_carlo.runs runs, monte_carlo.workers at a time
```

Each run is the full transient with inputs drawn from `monte_carlo.dispersions`: groups of
inputs, each with a relative or absolute error band, drawn uniformly. `"ALL_CDAS"` stands for
every CdA that carries flow, so an element added to a list is varied automatically. A group is
`"each"` (every input drawn separately, e.g. the CdAs) or `"shared"` (one draw moves the whole
group, e.g. the c* efficiency curve). Current bands are in `ASSUMPTIONS.txt`. The seed makes
the draws repeatable and run 0 is nominal. With `"thrust_control": "thrust"` the thrust targets
are re-solved with each run's engine (the controller holds thrust); `"flow"` flies nominal
flows. Results in `Outputs/MonteCarlo/`: every run's inputs and summary metrics
(`monte_carlo_runs.csv`), the statistics of each metric, i.e. the error bars
(`monte_carlo_stats.csv`: mean, std, 5/50/95 %), and time-history envelopes
(`monte_carlo_envelopes.png`).

## Plots

`Backend/Hopper_Plots.py` makes every graph in the old sim's output deck
(`outputs needed.pptx`) from the saved results, plus more. It can also plot any node, element,
tracked quantity or component attribute on request. `Hopper.py` makes the `outputs.plots` set
at the end of each run.

### The deck's graphs

| Deck slide | Old sim figure | Plot name |
|---|---|---|
| 2, 9 | Network Simulation: Nodes (pressure, temperature, mass, density, fill level) | `nodes` |
| 3, 11 | Ox Feed System pressure ladder (mdot, dP, CdA, State) | `ox_ladder` |
| 4, 11 | Fuel Feed System pressure ladder | `fuel_ladder` |
| 5 | Pressurization System pressure ladder | `pressurant_ladder` |
| 6, 10 | Engine Transient Performance (Pc, thrust vs target, MR) | `engine` |
| 6, 10 | Injector Stiffness | `stiffness` |
| 7 | Nozzle Exit Pressure vs Throttle | `exit_pressure` |

`plot deck` makes all seven.

The other named plots are:
- `feed_nodes`: the node grid for every station along both feed lines;
- `ox_stations` and `fuel_stations`: station pressures, plus the OTV or FTV CdA as a
  percentage of fully open;
- `ox_pressure_stations`, `fuel_pressure_stations`, `lox_pressurant_stations`,
  `fuel_pressurant_stations`: pressure against component along each path at chosen times
  (`outputs.station_times_s`, `--times`; `null` = peak thrust);
- `regen`: fuel temperature rise and heat pickup; `copv_temperatures`: bottles against their limit;
- `pressurization`: COPV and tank pressures, and the bang-bang states;
- `tank_states`: gas temperatures and remaining masses;
- `startup`: the main-valve opening over the first 2 s.

### On request

```bash
python Backend/Hopper_Plots.py list                                     # every plot, node, element, track and component
python Backend/Hopper_Plots.py plot deck                                # the deck's seven graphs
python Backend/Hopper_Plots.py plot ox_ladder engine --profile-time     # named plots, timed like the old sim
python Backend/Hopper_Plots.py plot ox_pressure_stations --times 1 5    # pressure along the LOX feed at 1 s and 5 s
python Backend/Hopper_Plots.py nodes "LOX Tank" COPV --csv              # node grid for any nodes or node groups
python Backend/Hopper_Plots.py ladder OTV OINJ "Regen Jacket" --tmin 5 --tmax 8   # ladder for any elements, zoomed
python Backend/Hopper_Plots.py track "Chamber Pressure [psia]" "Thrust [N]"       # any tracked quantity
python Backend/Hopper_Plots.py component OTV mass_flow upstream_pressure          # any stored component attribute
python Backend/Hopper_Plots.py plot exit_pressure --mr 1.25             # exit-pressure sweep at another MR
```

Options for every command:

| Option | Effect |
|---|---|
| `--file` | Results file to read (default `Outputs/Hopper_Prop_System.h5`). |
| `--out` | Output folder (default `Outputs`). |
| `--units E` / `--units SI` | psi and deg F, as the old sim (default), or Pa and K. |
| `--csv` | Also write each plot's numbers to `<plot>.csv`. |
| `--profile-time` | Time from the start of the thrust profile, the old sim's clock. This model's clock starts 0.5 s earlier, with the main-valve ramp. |
| `--tmin`, `--tmax` | Time window. |
| `--name` | File name for a custom plot. |
| `--show` | Also open the plot windows. |

Node groups are `tanks`, `ox_feed`, `fuel_feed` and `all`. Element groups are `ox`, `fuel` and
`pressurant`. Bad names stop with a one-line error that lists the valid ones. From Python:

```python
from Backend.Hopper_Plots import HopperRun, make_plots
make_plots(HopperRun("Outputs/Hopper_Prop_System.h5"), ["nodes", "ox_ladder"], "Outputs", csv_out=True)
```

## Tests

```bash
python -m unittest discover -s tests
python tests/run_feature_validation.py     # short coupled run of the header, valve-travel and regen-curve options
```

## Numerics

`simulation.nonlinear_max_evaluations` limits each SciPy nonlinear attempt, so a hard step
retries at a smaller time step instead of spending thousands of iterations on one attempt. The
residual tolerance is unchanged; a step is never accepted because its budget ran out. See
`DESIGN.md` and `HANDOFF.md` for the design record and limitations.

## Also included a comparison folder to sanity check with old sims

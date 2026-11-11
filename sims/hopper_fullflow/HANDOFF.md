# Hopper FullFlow handoff — November 5, 9 and 10, 2026

This document explains the work delivered during this development session,
the resulting code, the evidence behind the changes, and the remaining work.
The original November 5 sections describe the tree published in **`996f8509`**;
the November 9 update describes the independent torch purge implementation, and the
November 10 update (first below) the flight sequences, torch coupling and purge regulator.
This handoff does not replace each run's saved inputs or hardware calibration.
Dates and human-readable commit times below use America/New_York (EST).
Timestamped result-directory names use UTC.

> **Later change (November 7, 2026):** the element-list work from October 27 to
> November 3 was merged back into this structure, including the removed features
> listed in section 2.3. Modules moved to `Backend/`, some inputs and defaults changed,
> and links to flat module paths below are out of date. See the first section of
> [DESIGN.md](DESIGN.md) and [README.md](README.md) for the current state.

Start with [README.md](README.md) for operating instructions,
[DESIGN.md](DESIGN.md) for the ongoing design record, and
[hopper_config.json](hopper_config.json) for the actual inputs.
The [complete commit file manifest](docs/handoff_commit_files.txt) lists every
changed file, including generated outputs, for all three delivered commits.

## November 10, 2026 update: flight sequences, torch in the flight sim, purge regulator

This update describes the uncommitted changes on top of `083d88d3` (Vedanth's
November 9 torch purge commit), made by the user (nayangandra) with Claude. All
paths are in `sims/hopper_fullflow`. 54 tests pass; the 30 s `52m_free` run solves
cleanly in about 15–20 min (1072 steps), and its results are in `Outputs/`.

### Start-up and shutdown sequences (flight sim)

The single "both main valves open at 0.3 s over 0.2 s" sequence is replaced by the
proposed start-up and shutdown sequences.

- `hopper_config.json` `sequence`: `main_valve_start_s` and `demand_ramp_s` are
  removed. New fields: `main_fuel_open_s` (0.30), `main_ox_open_s` (0.45),
  `main_ignition_s` (0.70, profile time 0), `ox_close_after_shutdown_s` (0.0),
  `fuel_close_after_shutdown_s` (0.1) and `main_valve_close_s` (0.2, GUESS).
  `HOPPER_DATA.py` names: `MAIN_FUEL_OPEN_TIME`, `MAIN_OX_OPEN_TIME`,
  `PROFILE_START_TIME`, `OX_CLOSE_DELAY`, `FUEL_CLOSE_DELAY`, `MAIN_VALVE_CLOSE_TIME`
  (the old `MAIN_VALVE_OPEN_TIME` / `MAIN_VALVE_RAMP_TIME` are gone).
- `Backend/Hopper_Prop_System.py`: `burn_timeline` returns the ignition time,
  the shutdown request (end of the profile) and the time both valves are shut.
  `target_schedule` ramps each flow target from its own valve's start to the first
  profile point at ignition (fuel lead), and closes OTV first, FTV 0.1 s later
  (fuel-rich shutdown). New helper `valve_window`.
- Callers updated: `Backend/Hopper_Plots.py` (`PROFILE_START`),
  `Backend/Hopper_Report.py` (run metadata), `isolated_sims/Hopper_Purge.py`
  (`purge_start_time` = shutdown request), `Old_sims_comparison/compare_models.py`,
  `tests/test_port.py`.
- Fuel-only flow (T+0.30 to T+0.45) burns at the chamber map's lowest MR (0.5);
  propellant still lights on arrival.

### Purge sequence schema

- `purge.sequence` entries are now `{branch, open_after_shutdown_s, duration_s,
  note}`; `start_after_shutdown_s` and `gap_s` are removed. Purges may overlap.
  Current inputs: OPV opens with the OTV close and FPV with the FTV close (T+0 and
  T+0.1 after the shutdown request), 2 s each. `PURGE_SEQUENCE` is a list of
  `(branch, open_after_shutdown, duration)`; `PURGE_START_AFTER_SHUTDOWN` and
  `PURGE_GAP` are gone. `purge_windows` / `schedule_purges` take the shutdown
  request time.

### Torch network solved in the flight sim

Vedanth's `torch.network` is now built into the flight transient
(`torch.enabled` = true). The "transient builder is not implemented" error is
removed; with `enabled` false the network is listed as not modelled, as before.

- New `torch.sequence`: `oxidizer_open_s` 0.08, `fuel_open_s` 0.19,
  `close_after_ignition_s` 0.055, `purge_open_after_ignition_s` 0.03,
  `purge_close_after_main_purges_s` 0.1. Each network valve's
  `open_time_s` / `close_time_s` is its stroke (the schedule ramp).
  `Backend/Hopper_Config.py`: `_torch_legs` builds element objects from the
  network (OTHX `HeatedLine` becomes a tube with its friction-equivalent CdA, no
  heating), `_torch_sequence` validates the timing. New data names `TORCH_LEGS`,
  `TORCH_SEQUENCE`. The purge leg's sink (`FIORF Outlet`) is handled by putting the
  fuel-leg elements after the tee (the torch fuel injector) in series with it.
- `Backend/Hopper_Prop_System.py`: `torch_windows` / `schedule_torch` set the valve
  schedules; `build_torch` builds the LOX and fuel legs from their feed junctions,
  the TFPV purge leg from the regulator outlet, a zero-volume torch chamber on the
  main chamber map, and the torch throat into the main chamber. The torch LOX and
  fuel count in the main chamber's flow and MR. The TFPV nitrogen pressurizes the
  torch and leaves through its throat, but is not added to the main chamber (the
  map would treat it as hot combustion gas; counting it made thrust about 1.8 %
  high).
- `Backend/Hopper_Components.py`: `torch_throat_flow`, ideal subsonic/choked
  throat flow in either direction, for the mass-weighted ideal-gas mix of torch
  products and purge nitrogen, linear inside +/- the nozzle band at zero dP, with
  pressures and temperatures clamped positive (a negative-pressure Newton trial
  stalled the solver after shutdown). The torch map input is clamped to the map's
  lowest pressure for the same reason.
- `Backend/Hopper_Hardware.py` rows and `tests/test_torch_config.py` /
  `tests/test_config.py` updated for an enabled torch.
- Result: the torch reaches 510–550 psia before main ignition with the current
  CdAs (the proposed sequence chart shows about 150 psia; the 0.0065 in² throat
  is the main restriction), and TFPV keeps the torch 17–27 psi above MCC Pc
  through the burn.

### Purge regulator: Aqua Environment 1026, shared capacity

- `PR-REG` is now the Aqua Environment 1026 dome-loaded reducing regulator:
  Cv 0.05 (0.849 mm², was 5 mm² GUESS), set 750 psig, outlet rise 17 psi per
  1000 psi inlet drop (datasheet). New `Regulator` inputs
  (`Backend/Hopper_Elements.py`, `Backend/Hopper_Config_Fields.py`):
  `outlet_rise_psi_per_1000_psi_inlet_drop`, `reference_inlet_pressure_psia`
  (null = the source's starting pressure) and `full_open_droop_psi` (10 psi
  GUESS, not on the datasheet).
- `Backend/Hopper_Prop_System.py` `build_regulators`: with a CdA, the regulator is
  a proportional valve (shut at its target, fully open `droop` below it) feeding
  one outlet node that every purge and relief on it shares, so their combined
  demand is capped by its Cv. The outlet node stores the gas in the purge lines up
  to their solenoids (`regulator_outlet_volume`); its gas is throttled at constant
  enthalpy. Without a CdA it stays ideal (`min(source, set)`). `build_branches`
  takes `finish=False` so the torch purge can draw from the outlet before
  `finish_regulators` closes its mass balance. An exact-set regulator on a
  zero-volume outlet stalled the solver when a purge opened or closed on a fully
  open regulator; the droop model and outlet volume fixed it.
- Result: the 1026 is fully open during the shutdown purges, its outlet sags from
  about 812 to about 270 psia, and the fuel purge reaches 82 % of its requirement.

### Sizing tools

- `isolated_sims/Hopper_Purge_Sizing.py`: the blowdown steps overlapping purges
  together, and they, plus the torch purge (open throughout, to ambient), share the
  regulator capacity (`_solve_outlet`). `simulate` takes `regulator_cda` and
  `required_mass_flow` overrides. CLI `--order` is replaced by `--delays`. The
  "no minimum CdA" message now points at regulator capacity.
- New `isolated_sims/Hopper_Regulator_Sizing.py`: minimum regulator CdA / Cv versus
  set pressure for the purge requirement (or `--mass-flow` g/s per purge), with
  the path limit (no regulator is big enough) and the COPV dropout limit. Writes
  `regulator_sizing.csv`, `regulator_sizing_sequence.csv` and
  `regulator_sizing.png`; called from the report. Current result: sets from about
  748 to 1156 psig work with Cv ≥ 0.17 (2.85 mm²); the 1026 is too small at every
  set.

### New isolated tools and plots

- `isolated_sims/Hopper_Hotfire.py` with `isolated_sims/hotfire_config.json`:
  test-stand hotfire with constant-pressure run tanks (1000 psia), the vehicle
  engine, OTV/FTV, injector and regen, the torch legs into a zero-volume torch
  chamber and throat, and the purge regulator, relief and purges from a constant
  GN2 supply (2000 psia GUESS). Default profile: 100 → 500 lbf ramp over 2.5 s.
  Writes `Outputs/Hotfire/<UTC>/` (history CSV, per-graph profiles, maxima,
  summary and plots; the raw `.h5` is git-ignored). Its torch purge does not add
  gas to the torch, unlike the flight sim. Inputs documented in
  `docs/hotfire_inputs.md`.
- New `isolated_sims/Hopper_Results_Page.py`: builds `Outputs/results_page/`
  (`flight_overview.png`, `free flight simulation.csv`, `flight_maxima.csv` and
  `results_page.md`, the Notion results page with every input table) from the
  latest run.
- `Backend/Hopper_Plots.py`: new named plots `startup_sequence` and
  `shutdown_sequence` (torch and MCC Pc with every sequenced valve's opening).
- `ASSUMPTIONS.txt` and `.gitignore` (`Outputs/Hotfire/**/*.h5`) updated.

### Removed on November 11

The isolated sizing calculators were deleted (judged not good enough to keep):
`Hopper_Purge.py` (purge check), `Hopper_Purge_Sizing.py`, `Hopper_Regulator_Sizing.py`,
`Hopper_Vent.py`, `Hopper_Drain.py` and `Hopper_Tank_Pressure_Sweep.py`, with their report
hooks, their config inputs (`vent`, `lox_drain`, `tank_sweep`, and `purge`'s
`required_fraction_of_peak_flow`, `back_pressure_psia` and `sizing`), their `HOPPER_DATA`
names, their parameter-workbook rows and their outputs. The descriptions of them above are
historical. The run report now writes only `run_summary.csv` and `cda_table.csv`. Kept:
`Hopper_Hotfire.py`, `Hopper_Results_Page.py` and Vedanth's `Hopper_Torch_Purge.py`.

### Open items

- Pick a purge regulator with Cv ≥ ~0.2, set between about 760 and 800 psig (under
  the ORV-1 825 psig crack), or stagger the purges.
- Torch CdAs or throat versus the intended ~150 psia torch Pc.
- COPV gas reaches about 210 K, below the 220 K GUESS limit; that is the adiabatic
  gas temperature, not the wall. A wall heat-transfer model would answer it.
- Not modelled: spark plug, the 100 ms torch Pc check, OTHX heating, regulator dome
  loading; the chamber map is used before both propellants arrive.
- Hotfire line lengths are copied from the vehicle.

## November 9, 2026 update: torch topology and independent purge model

This update describes the changes after remote `4797ba36` (November 9 main-purge
work), building on `b3fe8cc4` (November 7 restoration). The November 5 sections
below remain historical. Current module paths and behavior are in README/DESIGN.
The publication commit contains the local changes in `sims/hopper_fullflow` only;
the separate NESS checkout is outside this work.

### Delivered implementation

- `hopper_config.json` declares the complete torch topology: FTV Inlet -> FIMV ->
  FICK -> FIORF -> fuel injector -> torch chamber -> throat -> main chamber;
  OTV Inlet -> OIMV -> OICK -> OIORF -> OTHX -> oxidizer injector -> the same
  chamber. Regulated GN2 flows through TFPV (renamed from FIPV), FIPCK and FIPORF
  into **FIORF Outlet**, downstream of FIORF and upstream of the shared injector.
- `Backend/Hopper_Torch_Config.py` validates topology, units/geometry, unknown or
  duplicate connections, command schedules, and OTHX's prescribed dT-versus-time.
  OTHX time zero is torch opening; linear interpolation holds endpoint values.
  This is input/helper support, not a solved hot oxidizer heat exchanger.
- `Backend/Hopper_Config.py` and `HOPPER_DATA.py` expose the normalized topology;
  hardware inventory and parameter-workbook descriptions include all torch parts.
  Legacy `torch.branches` remains replayable. Enabling the explicit network in
  the main vehicle solver raises an error instead of silently using the old
  direct-main-chamber approximation. Keep `torch.enabled=false`.
- `Backend/Hopper_Torch_Purge.py` independently solves common compressible gas
  flow and individual station pressures through the purge restrictions and
  shared fuel injector. After pressure release, BDF integrates log chamber mass
  and internal energy plus cumulative inlet, outward-throat and reverse-flow
  masses. Positive throat flow is torch-to-main; reverse flow is permitted.
  Geometric throat area is used without an extra discharge coefficient. The
  near-zero pressure-drop regularization is an explicit numeric input.
- `torch.standalone_purge.sequence` owns the complete event timeline. Each event
  can change valve commands and main pressure and release imposed torch pressure.
  Actual valve travel follows separate full-stroke open/close times; partial
  reversals are supported. Commands hold between events. Legacy per-component
  command archives still replay, but mixing both formats is rejected.
- `isolated_sims/Hopper_Torch_Purge.py` is the independent CLI, with `--config`,
  `--out`, `--check-inputs`, `--diagram-only` and `--no-plots`. It archives exact
  input JSON, replayable config/assets, derived settings and the resolved event
  timeline. Histories, component station pressures, steady-equilibrium checks
  and summaries are exported alongside PNGs.
- `Backend/Hopper_Torch_Diagram.py` generates SVG/PNG topology drawings and a
  component-edge CSV. The shared tee is explicit; dashed propellant legs are
  configured only and the solid purge path is solved. Dimensions are imperial,
  CdAs are mm2, and throat area/diameter are in2/in.
- Purge-only pressure ladders span opening, steady flow, closing and the tail.
  PNGs include a torch-to-main ladder zoom, station histories, pressure-margin
  history, startup transition, flow/valve history and component diagram. CSVs
  supplement the PNGs; they do not replace them.

### Current inputs versus archived comparison inputs

The current checked-in config preserves the user's later **fuel injector CdA
0.8051 mm2** and **FIPORF CdA 50 mm2**. FIORF/OIORF are provisionally 8 mm2.
FIMV/TFPV travel is 0.08 s and OIMV travel is 0.10 s. OTHX is 60 in long,
0.25 in OD with 0.065 in wall (derived ID 0.120 in), with its user dT curve.
Torch volume is 0.0019 L; throat is 0.006504 in2 (about 0.0910 in diameter).
Initial/imposed torch pressure is 150 psia; initial/source/main-reservoir
nitrogen temperatures are assumed 293 K with gamma 1.4. The constant source
resolves to 764.8 psia from the initial regulator/bottle configuration.

The base config still uses the example 5-second purge: open command at 2 s,
close at 7 s, finish at 7.5 s. At 2 s the feeds receive close commands, main
pressure steps to 300 psia and imposed torch pressure is released. This is not
30 seconds: the requested 30-second cases were separate archived configs.

The archived later-input run `Outputs/Torch_Purge/20261110T005911_735990Z/`
(UTC directory; November 9 EST) used the 0.8051 mm2 injector. Its summary reports
9.5020 g/s steady purge, 315.6931 psia torch at 300 psia main, and 47.9363 g
consumed for the example command window including travel/tail. This is distinct
from the smaller-injector results discussed earlier in the session.

The two requested 30-second cases live under
`Outputs/Torch_Purge/comparison_30s_20261109T183407_098289Z/`. Their saved inputs
use **0.1051 mm2 fuel injector CdA** and 50 mm2 FIPORF. TFPV commands open at
2 s and close at 32 s; integration ends at 32.5 s, including 80 ms closing travel.

- Main 300 psia: 1.2860077 g/s; total source gas 38.6758915 g; steady torch
  300.2863972 psia, or +0.2863972 psi.
- Main 160 psia: 1.2860077 g/s; total source gas 38.6759284 g; steady torch
  160.5370870 psia, or +0.5370870 psi.
- Ideal instantaneous 30-second consumption is 38.5802304 g in either case.
  The actual command window plus closing travel produces the slightly larger
  integrated totals. Maximum mass-balance error was below 3.3e-11 kg.
- Both cases choke at the fuel injector with fixed upstream pressure and
  temperature; changing main pressure therefore does not change steady mdot.
  The torch throat itself is unchoked. A comparison PNG, per-case PNGs, input
  archives, numeric histories and `comparison_summary.json` are included.

Earlier timestamped directories retain intermediate configurations, so do not
assume they match today's base config. The first two input-only directories are
failed development attempts, not completed simulations. Always read each run's
`input_config.json`, `run_config.json` and `run_summary.json` together.

### Running and validation

From this folder, with dependencies installed:

```bash
python isolated_sims/Hopper_Torch_Purge.py
python isolated_sims/Hopper_Torch_Purge.py --check-inputs
python isolated_sims/Hopper_Torch_Purge.py --diagram-only
python isolated_sims/Hopper_Torch_Purge.py --config <saved-run>/run_config.json
python -m unittest discover -s tests
```

On this workstation, `./.venv-validation/Scripts/python.exe` is the working
interpreter. Automated checks cover FullFlow agreement for choked/subcritical
orifice flow, equal restriction flows, closed/reverse supply paths, full/partial
valve travel, event release/backpressure, mass conservation, timestep refinement,
central timeline validation and archive compatibility, steady throat balance,
throat-area sensitivity, topology tee edges and OTHX interpolation. The cold-gas
regression fixture pins its physical inputs rather than relying on the editable
fuel-injector CdA. The complete suite contains 54 tests. PNGs were rendered and
visually checked; the two comparison cases were re-exported from their saved
histories without re-solving or changing their physical inputs.

### Physical limits and agreed interpretation

This is **cold ideal-gas purge only**, not torch combustion or a coupled vehicle
purge solve. No residual propellant flow through closing feeds, liquid clearing,
phase change, OTHX heating, wall thermal dynamics, hot-main-gas chemistry, trapped
line storage, check cracking pressure, regulator droop or COPV blowdown is modeled.
Even main-to-torch backflow is represented as the configured cold nitrogen gas.

The 150 psia pre-release phase is imposed; it does not validate torch firing
performance. Source consumption is nitrogen admitted through the purge inlet,
not gas entering from the main chamber. The steady pressure margin is an output,
not an independently imposed pressure or a demonstrated ingestion protection limit.
A static pressure boundary does not model chamber oscillations, local turbulence,
diffusion, simultaneous inward/outward throat exchange or hot-gas penetration.

The session discussed +/-5% **peak** main-pressure amplitude (15 psi at 300 psia,
8 psi at 160 psia), with +/-1% and +/-10% sensitivity cases, as unvalidated
exploratory assumptions. No fluctuation waveform, RMS/peak conversion, frequency
spectrum, safety factor or purge acceptance criterion has been implemented.
Comparing a mean pressure margin with those amplitudes alone does not establish
reverse-flow mass or thermal protection. The proposed 3x supply-pressure scaling
was an estimate, not a simulated/qualified case. No hardware topology bypass or
injector redesign was implemented.

Next work: characterize local main pressure versus time at the torch outlet;
define allowable reverse-flow mass/thermal exposure and minimum outward flushing;
then evaluate a transient backpressure envelope. Add gas/thermal realism or full
combustion/OTHX coupling only as required for those acceptance criteria. There is
currently no validated minimum purge mdot or automatic hardware pass/fail.

## Contents

1. [Objective and final state](#1-objective-and-final-state)
2. [Commit history and publication](#2-commit-history-and-publication)
3. [Architecture and execution](#3-architecture-and-execution)
4. [Changes in propulsion modeling](#4-changes-in-propulsion-modeling)
5. [JSON configuration migration](#5-json-configuration-migration)
6. [Plots, diagram, reports, and run archives](#6-plots-diagram-reports-and-run-archives)
7. [Runtime investigation and fixes](#7-runtime-investigation-and-fixes)
8. [Final checked-in inputs](#8-final-checked-in-inputs)
9. [File-by-file implementation guide](#9-file-by-file-implementation-guide)
10. [Validation and its limits](#10-validation-and-its-limits)
11. [Running, replaying, and extending the code](#11-running-replaying-and-extending-the-code)
12. [Known limitations and follow-up priorities](#12-known-limitations-and-follow-up-priorities)

## 1. Objective and final state

The user's objective is a configurable propulsion-system simulator: supply a
thrust history and hardware/initial-condition inputs, then predict pressures,
temperatures, flows, inventories, and engine performance over time. Component
choices and geometry should be easy to change and visible in one input file.

The delivered implementation adds:

- Pressure versus named component station at selected times.
- Prescribed fuel-regen temperature rise versus chamber pressure, including
  downstream liquid-property updates.
- A generated component/source schematic.
- Named thrust profiles and a parameterized ramp; sustained demonstrated-thrust
  reporting with explicit interpretation of capacity evidence.
- Separate valve opening and closing travel, independent controller sampling,
  pressure deadbands, and retry-safe control state.
- Individual COPV inventories/thermal states and minimum-temperature flags.
- Tank propellant residual reports.
- Engine c-star and thrust-coefficient efficiency functions or CSV tables.
- A structured JSON input file, unit conversion/validation, and per-run archives.
- Cleaner plot panels, imperial diagram geometry, a Git Bash launcher, and a
  bounded nonlinear solve attempt to recover from an early timestep stall.
- Unit tests, coupled feature checks, and design/user documentation.

The model remains a single-phase feed system with a quasi-steady, zero-volume
chamber. This work did not add a full combustion transient, ignition model,
two-phase tank depletion, or a wall-resolved regenerative-cooling model.

## 2. Commit history and publication

### 2.1 `31d3aad8` — initial propulsion features

- Full hash: `31d3aad81e48a9d8b0d32b3e5a20529d6382fba1`.
- Time: November 5, 2026, 03:17:54 EST.
- Subject: **Add configurable Hopper propulsion dynamics and reporting**.
- Parent: `748458d4fd8d242603099c8ef4dbecb258820654`.
- Git totals: 50 files changed; 6,508 insertions and 763 deletions.

This commit delivered the initial station plots, prescribed regen heating,
box diagram, ramp/profile selection, valve travel and sampled tank control,
individual COPVs, residual and maximum-thrust reports, efficiency inputs,
tests, and design document. Inputs still lived in `HOPPER_DATA.py` at this stage.
It also committed sample reports/plots and ignored local environments,
validation scratch files, Python caches, and raw propulsion HDF5 histories.

### 2.2 `9b85078b` — JSON, archives, presentation, and runtime recovery

- Full hash: `9b85078b4dfe0477fbd2c2c84cf415641a790d8c`.
- Time: November 5, 2026, 05:27:41 EST.
- Subject: **Move Hopper inputs to JSON and improve run diagnostics and solver retries**.
- Parent: `31d3aad8`.
- Git totals: 143 files changed; 29,203 insertions and 3,883 deletions.

This commit replaced editable Python vehicle inputs with JSON, added the loader
and field bindings, archived each driver's inputs/assets, redesigned the
schematic, removed Hdot/quality plot panels, and added the Git Bash launcher.
It includes the fuel-pressure correction, progress/total-runtime reporting,
the nonlinear-attempt budget, additional tests, two trajectory CSVs, and the
three committed timestamped runs described in section 10.

Some parameter changes were user edits included in the requested "commit all"
operation, rather than calibrated defaults selected by the implementation.
For example, the final ramp is 2.5 seconds, the pressurant orifices are roughly
6 mm², and throttle maximum CdAs are 50 mm². Do not confuse those final settings
with the earlier 10-second performance benchmark.

### 2.3 `996f8509` — FullFlow-only replacement on remote `main`

- Full hash: `996f850999e15ed25b5a0cc16f338024f860b419`.
- Time: November 5, 2026, 14:02:57 EST.
- Subject: **Replace FullFlow with local JSON-configured simulator**.
- Parent: remote-main commit `f23441ad2e8aeb06ef2752e80d2eb5d2dfd499d1`.
- Git totals: 225 files changed; 28,468 insertions and 5,788 deletions.

The publication history matters:

1. Local `main` was two commits ahead of, and 212 commits behind, remote `main`.
2. Remote FullFlow had independently been reorganized into `Backend/`,
   `Outputs/`, `thrust_curves/`, and `isolated_sims/`, with additional features.
3. A merge was tested in an isolated worktree and produced 23 conflict paths.
   The merge was aborted. It was never published.
4. At the user's initial selection, the two local commits were published to
   `hopper-json-config` for review.
5. The user then explicitly instructed that remote FullFlow be fully replaced
   by the local version, with no local changes outside FullFlow pushed.
6. A new commit was built on current remote `main`, replacing only
   `sims/hopper_fullflow` with the exact tracked tree from `9b85078b`.
7. The commit was pushed normally to `main`; no force push or remote-history
   rewrite was used. Local `main` was aligned to that commit afterward.

Consequently, **`31d3aad8` and `9b85078b` are not ancestors of the final `main`
commit**. Their resulting FullFlow tree was copied into `996f8509`, whose parent
is the previous remote main. Git verified that the FullFlow subtree hashes were
identical:

```text
9b85078b:sims/hopper_fullflow = 0f1f62f3bb38fd3a95c122e098f1ae30db28c658
996f8509:sims/hopper_fullflow = 0f1f62f3bb38fd3a95c122e098f1ae30db28c658
```

Everything outside FullFlow in `996f8509` matches its remote parent. Unrelated
local changes in the Simulink file, general-fluid-network settings/result, and
3-DOF test script were discarded as requested. The separate untracked `NESS/`
folder, containing a nested repository/environment, was explicitly kept.

This was an **intentional folder replacement**, not an integration of the two
FullFlow implementations. The remote-only Monte Carlo driver, spreadsheet
parameter machinery, alternate backend/report/schematic modules, isolated
drain/purge/vent/pressure-sweep scripts, and their output organization were
removed from the final FullFlow tree. Their earlier versions remain in Git
history. Restoring those features requires a deliberate port into the JSON-based
implementation; they are not delivered features of this handoff.

The manifest's deletion/rename entries make that replacement explicit.
It does not imply deletion of unrelated remote repository features.

## 3. Architecture and execution

### 3.1 Config-to-results pipeline

```text
hopper_config.json or --config / HOPPER_CONFIG
    -> Hopper_Config.load_config + Hopper_Config_Fields conversions
    -> HOPPER_DATA compatibility namespace (SI values and hardware objects)
    -> Hopper_Run_Config.prepare_run (driver solves only)
    -> optional chamber-map generation
    -> steady thrust-target solves
    -> hardware diagram and target CSV
    -> transient network construction + liquid property tables
    -> bounded-attempt FullFlow transient solve
    -> HDF5 metadata/config
    -> summary, residual CSV, performance CSV, plots, timing
```

`Hopper.py` is a top-level executable script. Importing it executes argument
parsing and work; it is not a library-safe `main()` function yet.
Library callers should use the targeter/network/report functions directly.

Several modules use `from HOPPER_DATA import *`, and some function default
arguments capture values at import time. The driver deliberately sets config
selection and rewrites paths to archived assets **before importing the physics
modules**. Avoid changing configuration globals after those imports and
expecting every consumer/default argument to update automatically.

### 3.2 Network topology

Each individual COPV supplies a common gas header through a configurable
connection restriction. The header supplies oxidizer and fuel pressurant legs.
Each leg connects to its tank's ullage. The liquid paths connect tanks to the
common chamber/nozzle.

Storage volumes integrate mass and, for gases, energy. Liquid feed junctions
and the chamber are zero-volume pressure nodes. Liquid tubes include inertia
and friction. Restrictions may include passage inertia when a length is given.
Gas lines can use a lumped equivalent restriction or a dynamic FlowTube.
Consecutive gas restrictions are grouped into a SeriesCdA/CompressibleOrifice.

All state histories are coupled by the implicit transient solver; a configured
thrust history requests flows, but the hardware can saturate and underdeliver.

### 3.3 Engine interpretation

The equilibrium chamber map contains ideal gas-product temperature, gamma,
gas constant, and density as functions of chamber pressure and mixture ratio.
It is generated with tank-temperature reactants. It does **not** directly store
the empirical efficiency curves introduced here.

`EnginePerformance` combines that gas map with the nozzle and separate
efficiency lookups:

```text
real chamber mass outflow = ideal nozzle mass flow / eta_cstar
delivered thrust = eta_cf * (ideal momentum thrust + ideal pressure thrust)
```

The targeter solves required pressure at the configured target O/F for each
profile point, then splits total flow into fuel and oxidizer. The transient
uses the same engine block and computes actual mixture ratio from chamber
inflows, with a small flow floor and map-range clipping for numerical stability.
The thrust-capacity report additionally checks unclipped MR against map bounds.
Startup nozzle flow has a small-pressure-drop regularization near ambient.

There is no finite chamber gas inventory, ignition delay, combustion instability,
or wall/combustion thermal state in this model.

## 4. Changes in propulsion modeling

### 4.1 Thrust profiles, sequencing, and demonstrated capacity

The initial feature work added validated ramp generation and profile validation.
Times must start at zero and be unique; thrusts must be finite and positive.
Both increasing and decreasing positive-thrust ramps are supported.

The JSON form selects `ramp` or a named inline profile. The ramp has start/end
thrust in lbf, duration in seconds, and a number of sampled points. Named
profiles use `[time_s, thrust_N]` pairs. The original 29.7-second `1dof` profile
remains available.

`linear` interpolation interpolates the requested flow targets between solved
points. `step` holds demands with the configured step-smoothing interval.
Neither option replaces physical valve travel. Startup/shutdown demand ramps
and post-shutdown tail-off are separate sequence inputs.

The maximum-thrust report is deliberately run-specific. It finds the largest
minimum thrust maintained throughout a full configured hold interval, using
saved samples and excluding invalid inventory/Pc/MR samples. It therefore
rejects a short spike as a sustained result. A separate test looks for sustained
valve saturation plus target shortfall. Without that evidence, the achieved
thrust is reported as a lower bound, not a proven system maximum.

This does not implement an independent global capacity search or steady
optimization over every tank state, mixture ratio, or hardware combination.
Save interval and hold time affect the sampled sustained-thrust calculation.

Two trajectory CSVs were added with the user's files:
`thrust_profile_10m_tethered.csv` and `thrust_profile_52m_free.csv`. Their columns
include reference/simulated thrust and altitude. **The JSON loader does not
automatically load thrust CSV files.** Choose the intended thrust column and
convert it into a named inline `time_s_thrust_n` profile, or add a validated
CSV-thrust input feature. Efficiency CSV support is a separate feature.

### 4.2 Physical valve travel and sampled bang-bang control

`ValveTravel` tracks normalized opening with separate full-stroke opening and
closing times. Partial travel takes proportional time. A zero travel time is
ideal instantaneous movement. Position is bounded to `[0, 1]`.

Travel is computed from the previous **accepted** state and elapsed time.
Repeated residual evaluations and rejected solver attempts do not repeatedly
advance a valve. This is essential for nonlinear solves and timestep retries.

`SampledBangBang` holds the last command between controller samples and applies
pressure hysteresis using accepted-state pressure. Sample times are solver
breakpoints. Its interval is independent of valve travel and requested solver
step. At a sample boundary, command evaluation is tied to the accepted state
at the start of the following integration interval.

Tank control opens below `setpoint - low_deadband`, closes above
`setpoint + high_deadband`, and holds its state within the band. Each tank has
an initial controller command. In the gas builder, physical solenoid position
starts closed and responds to that command through travel. Thus an initial
command of one can create a pre-burn pressurization pulse even at set pressure.

Throttle control now solves a requested effective area separately from physical
area. A rate-limited position determines actual CdA, limited by maximum CdA.
Actual hardware capacity can then prevent a target from being met. Scheduled,
bang-bang, throttle, and check valves expose open/close timing where supported.
Check valves also enforce forward delivery. Fixed restrictions do not acquire
actuator timings.

For throttle valves, the legacy ladder `State` display can still be the clipped
target-flow/full-scale-flow ratio. It is not necessarily measured actuator
position. Use the explicit `Opening` track or actual CdA/max-CdA to inspect
physical travel.

The closed limit of an optional inertial restriction was guarded to avoid
division by zero when effective CdA is zero. Gas series restrictions combine
each component's own opening with its CdA before calculating total effective
area, preserving partial-position effects.

### 4.3 Regen temperature rise and downstream properties

The original regen was a hydraulic orifice. `Regen` retains hydraulic CdA and
optional passage length, and adds a prescribed `delta_t_curve` evaluated at Pc.
In JSON this is `pc_psia_delta_t_k`.

```text
fuel outlet temperature = fuel inlet temperature + prescribed delta T(Pc)
prescribed heat pickup = mass flow * (outlet enthalpy - inlet enthalpy)
```

Pressure interpolation is linear; endpoints are held. The curve is a
**temperature rise in kelvin**, not an absolute outlet temperature. A decreasing
rise with Pc is accepted if that is the supplied curve. A first point at
150 psia also holds its rise below 150 psia, including unlit startup; there is
no automatic zero-Pc/zero-heating assumption.

The liquid feed builder updates downstream temperature and fluid properties
after each regen element, including injector inlet temperature/density. It
tracks inlet/outlet temperatures, rise, and heat pickup. Multiple downstream
elements use the updated context.

Direct liquid property calls are cached in memory. A temperature property table
is precomputed over the configured rise range, at the tank set pressure, and
interpolated inside nonlinear iterations. This avoids expensive enthalpy
integrations at every residual evaluation. The table currently rebuilds in a
new Python process; it is not a persistent disk cache.

The model does not solve coolant-wall heat transfer, wall inertia, heat flux,
reverse-flow heat transport, or feedback of the heated reactant into the
equilibrium chamber map. It assumes liquid validity over the configured range.

### 4.4 Tanks, inventories, and ullage collapse

Tank initialization uses **tank volume plus initial liquid mass**, not a target
ullage percentage. Initial liquid volume is mass divided by the selected liquid
density. Initialization checks that positive ullage remains.

Tank liquid temperature is fixed in the tank model. Ullage energy includes
convection to the liquid surface, with area `pi * radius^2` and the configured
HTC. Those inputs expose the existing ullage-collapse treatment. A separate
empirical mass-collapse multiplier was not added.

Residual reports distinguish initial, profile-end, and final liquid mass,
consumed mass, and percent of initial mass remaining. They refer to **tank
inventory**, excluding trapped feed-line inventory in zero-volume junctions.
Negative final inventory is flagged; nonpositive inventory samples are excluded
from demonstrated-thrust intervals. The solver does not switch into a physical
two-phase depletion model or automatically stop at empty tanks.

### 4.5 Individual COPVs and positive-volume header

The combined source was expanded into individually configured bottles with
their own initial pressure, temperature, volume, outlet CdA, mass, and energy.
Each bottle has an optional minimum-temperature operating limit.

The user explicitly selected **flag the temperature limit and keep solving**.
The report records minimum observed temperature, limit crossing, and first
crossing time. It does not clamp temperature or shut off pressurant delivery.

Bottle connections feed a common adiabatic header. A zero-volume header became
underdetermined when both delivery valves closed, so a positive gas-storage
volume was added. The current 0.05 L header and 100 mm² bottle connections are
model inputs that need measured plumbing data; they strongly affect stiffness.
Reverse delivery into the bottles is blocked.

Legacy aggregate COPV **pressure and temperature** tracks now describe the
header. Aggregate COPV **mass** sums individual bottle masses; header mass is
tracked separately. Downstream analysis must respect this distinction.

### 4.6 Engine efficiencies and the function-lookup correction

Both the targeter and transient engine use the same c-star/Cf efficiency
evaluation. The original Python-input stage allowed functions or measured CSVs.
The JSON stage represents functions declaratively as bounded linear functions
or Pc/MR polynomials, or selects a CSV file.

- Linear mode interpolates low to high over a pressure interval and holds edges.
- Polynomial coefficient `[i][j]` multiplies `Pc_psia^i * MR^j`; a single
  coefficient is a constant function. Values are validated during evaluation.
- CSV columns are `pc_psia,efficiency`, optionally with `mr`. A Pc/MR table must
  form a complete rectangular grid. Interpolation is linear/bilinear with held
  edges. Duplicate points, missing cells, nonfinite data, and efficiencies
  outside `(0, 1]` are rejected. Paths resolve relative to the config file.
- A non-null CSV can override a linear-mode function. JSON does not execute
  arbitrary Python function text or expressions.

The old low/high formula used `max(1, fraction)`, which held the high value at
low pressure and could rise without a cap. It was corrected to a fraction
clamped to `[0, 1]`. This is a real behavior correction; old comparisons cannot
be assumed numerically equivalent.

Another important fix is `FunctionLookup` in `Hopper_Control.py`. FullFlow 3 can
reuse a returned dictionary by calling `dict.update` with new inputs. That
updates input keys without recomputing calculated outputs, freezing functions
such as temperature rise or efficiency. The adapter sets `reuse_existing=False`
on the Lookup object, forcing function reevaluation. Passing it as an ordinary
constructor keyword would instead pass it to the wrapped function.

### 4.7 Geometry and fluid passage inertia

Liquid and gas lines accept OD and wall thickness and derive ID; alternatively
they accept ID alone. If all dimensions are specified, they must agree.
Invalid/nonpositive geometry and tube lengths are rejected. A tube edit does
not automatically resize an independent valve CdA.

Liquid tube inertia and optional restriction-passage inertia represent the
acceleration of fluid in a finite passage. They are different from actuator
travel. Gas lines use the configured inertia switch; their lumped CdA path
does not resolve a full distributed tube transient. Positive liquid-line
height change means the outlet is higher and gravity opposes forward flow.

## 5. JSON configuration migration

### 5.1 Input ownership and units

`hopper_config.json` is the editable user interface. `HOPPER_DATA.py` now loads
that file and exports SI-valued constants/hardware objects for existing modules.
`Hopper_Config_Fields.py` holds field-to-symbol bindings and fixed conversion
factors, not vehicle settings.

The loader checks schema version, required/unknown keys, duplicate JSON keys,
finite numbers, selected profiles, geometry, names, and key timing/range inputs.
Optional constructor defaults are materialized in the resolved config. This
is a custom validator, not an exhaustive external JSON Schema for every possible
type or physical operating condition.

Unit suffixes are authoritative: psia for absolute pressure, psi for pressure
differences, inches for geometry, in² for engine areas, mm² for CdA, kelvin for
temperature/rise, liters for volumes, kg for mass, and seconds for time.
The internal physics still uses SI. JSON `null` disables a temperature check
or represents omitted geometry as appropriate. JSON comments are unsupported.

### 5.2 Configuration sections

- `simulation`: requested maximum/save step, residual tolerance, Jacobian,
  per-attempt nonlinear budget, solve/rebuild switches, console output.
- `files`: chamber map, raw history filename, results root.
- `outputs`: selected plots, units, CSV/show/profile-time options, station times,
  schematic switch, nozzle-separation reference.
- `fluids`: oxidizer, fuel, pressurant, ambient pressure.
- `sequence`: initial demand start, demand ramp, step smoothing, tail-off.
- `numerics`: flow floor, nozzle regularization, property-table points,
  controller/valve step-resolution factors.
- `engine`: O/F, throat and exit areas, nominal pressure, efficiency definitions.
- `chamber_map`: Pc/MR grid limits and resolution.
- `pressurization`: controller interval and header volume.
- `maximum_thrust`: sustained duration and tracking-error criterion.
- `thrust_profile`: selected profile, interpolation, ramp parameters, named pairs.
- `bottles`: separate source objects and limits.
- `tanks`: geometry, initial inventory/temperatures/pressure, HTC, control bands.
- `flow_paths`: ordered liquid-feed and pressurant components.

The loader supports built-in serialized `Tank`, `COPVSource`, `Line`, `GasLine`,
`Orifice`, `Injector`, `Regen`, `Valve`, `BangBangValve`, `CheckValve`,
`ThrottleValve`, and `CavitatingVenturi` objects. A custom Python component still
needs a code/loader extension; it cannot be inserted into JSON as executable code.

### 5.3 Selection and import behavior

`python Hopper.py --config path/to/case.json` selects a case for the driver.
`HOPPER_CONFIG` selects a case for standalone tools too. With neither selection,
the default is the file beside `Hopper_Config.py`. Relative asset/results paths
resolve relative to the selected config's directory, rather than shell CWD.

Changing a CSV/function efficiency does not require regenerating the ideal gas
map. Changing propellants, reactant temperatures, or map grid should trigger a
rebuild. Existing map provenance is hashed/copied, but semantic compatibility
with all changed input fields is not automatically proved.

## 6. Plots, diagram, reports, and run archives

### 6.1 Station snapshots and plot cleanup

Four new plots show pressure against source/component names:

- `ox_pressure_stations`.
- `fuel_pressure_stations`.
- `lox_pressurant_stations`.
- `fuel_pressurant_stations`.

Station labels identify the source followed by component outlets. Times can be
selected in config or via `--times`; `null` selects peak delivered thrust.
Gas pressures inside a grouped restriction are allocated from effective CdAs,
not independently solved compressible volumes; the plot labels that approximation.
Plot CSV export was updated to handle string station labels.

`regen` and `copv_temperatures` add thermal/operating-limit plots, giving 19
named plots in total. The original engine, ladder, startup, stiffness,
pressurization, tank-state, nozzle-exit, and station-time plots remain available.
Hdot and quality panels were removed from standard node/ladder plots and their
plot CSVs. Ladders now have four panels. Node grids show the remaining five
quantities and hide the unused sixth slot.

Raw enthalpy/quality tracks were retained for compatibility. The regen heat-pickup
plot also remains: it is the heating enthalpy difference, not the removed
absolute enthalpy-transport Hdot panel. Thermal exports add available regen/node
temperature information while retaining older-history fallbacks where implemented.
History readers check saved path names against the currently loaded hardware;
use the archived config when plotting an old run.

### 6.2 Component diagram

The initial diagram was a wide box-and-arrow P&ID. The later design uses compact
serpentine oxidizer/fuel paths, rounded cards, separate branch colors, right-angle
connections, and simple tank/valve/line/restriction symbols.

It includes sources, header, tank inventories/setpoints, ordered hardware,
CdAs, timings, and chamber/nozzle geometry. Every line reports OD/wall status,
ID, and length in inches. Missing OD/wall are marked unspecified, not fabricated.
Throat area is entered in in²; throat diameter is derived and displayed in inches.
The diagram is saved as PNG plus scalable SVG. Symbols/layout are schematic and
do not represent actual physical location, routing distance, or certified P&ID
symbology.

### 6.3 Reports and metadata

`Hopper_Report.py` writes run-specific profile timing, integration step, map
limits, initial tank masses, throttle maximum areas, bottle limits, and path
names into `hopper_metadata` in HDF5. The final implementation also embeds
`hopper_config`.

Reports produce `run_summary.json`, `run_summary.txt`, and
`propellant_residuals.csv`; a summary is printed. They include sustained and
instantaneous peak thrust, capacity-evidence interpretation, inventories,
temperature-limit crossings, incomplete-profile status, and excluded samples.
Old histories without metadata cannot reliably supply the new summary and
must be rerun. `invalid_inventory` uses final tank mass, while nonpositive
inventory samples are separately excluded from capacity evaluation.

### 6.4 Driver-created run archives

Every new driver solve creates a unique UTC timestamp directory under the
configured results root before the expensive solve. It contains:

- `input_config.json`: original supplied JSON text.
- `run_config.json`: config with defaults materialized and portable paths into
  this run's local assets.
- `resolved_inputs.json`: config, sampled selected demand, derived throat
  diameter/expansion ratio, integration cap, header initialization metadata,
  and asset SHA-256 hashes.
- `inputs/chamber_map.h5`: copied or generated map actually used by this run.
- `inputs/cstar_efficiency.csv` / `cf_efficiency.csv` when selected.
- Raw history, target/performance CSVs, reports, selected plots/plot CSVs,
  diagram files, and `run_timing.json` after successful completion.

The physics modules are imported after the driver switches paths to those copied
assets, so later changes to original CSV/map files do not change that run's input
files. Hashes are refreshed after optional map generation. Failed attempts can
leave input directories without a completed history/report; a directory's
existence is not proof of successful completion.

These guarantees describe the **driver workflow**. Low-level library calls and
the in-process feature test do not invoke the complete run-archive workflow.
Archives record numerical inputs/assets, not a frozen copy of every Python
source or a dependency lockfile. Preserve the Git revision and environment
alongside important results.

### 6.5 Tracked and ignored generated files

All three commits include generated outputs as listed in the manifest. Root
`Hopper_Results/` files are historical sample outputs; the timestamped folders
have their own input snapshots and are better provenance sources.

`.gitignore` excludes `Hopper_Prop_System.h5` by basename (including within run
folders), local `venv/` and `.venv-validation/`, Python caches, and
`validation_outputs/`. Thus a fresh clone contains committed plots/configs and
maps, but not the raw transient histories needed to regenerate arbitrary plots.
Keep/copy those large histories separately or rerun the archived config.
`.venv/` is supported by the launcher but is not explicitly listed in this
folder's `.gitignore`; check repository ignore rules before staging a new one.

## 7. Runtime investigation and fixes

### 7.1 Why the first feature work increased runtime

The previous editable requested step was 50 ms, but physical/controller
resolution now caps actual integration at 5 ms for the checked settings.
That alone increases the step count by roughly ten times. Separate bottles
plus header add mass/energy states and nonlinear gas-property work. A small
header and large outlet areas create fast gas dynamics. A 101-point ramp also
requires more steady target solves than a shorter profile list.

The JSON parser, archive copies, and drawing/plot changes did not create the
early nonlinear timestep stall. They do contribute startup/output work, but
the measured bottlenecks were fluid-property construction and hard transient
solve attempts.

### 7.2 Input errors found during investigation

- Regen entries of 346/339 K were initially entered as rises, giving approximately
  639/632 K outlets with 293 K inlet. The user confirmed they had entered outlet
  temperatures by mistake, then changed the curve to rises of 51/44 K. Those
  final values mean 344/337 K outlets. The implementation retained the latest
  chosen curve rather than silently changing it to 53/46 K.
- Fuel tank initial pressure was 5,004 psia with a 500 psia setpoint and bottles
  initially at 4,503.2 psia. The user explicitly confirmed a typo; it was set to
  500 psia.
- A positive LOX-line height of 48 in was identified as a four-foot upward run,
  not gravity assistance. The final committed configuration has height change
  zero. That later parameter edit is distinct from the runtime benchmark.

The final snapshot also has user-selected geometry, orifices, throttle limits,
and shorter ramp settings; these are recorded in section 8 and run configs.

### 7.3 Heated-property construction delay

Profiling the valid warm-fuel case showed about 59 seconds of property-table
construction under profiler overhead, mostly enthalpy integrations in
ThermoProp/RocketProps. An unprofiled checked build took about 26 seconds;
other runs varied. These are setup costs before timestep integration.

Added console messages identify network/property construction, show table
temperature range/count/pressure, print periodic table progress, and announce
when the transient actually starts. The driver prints initial tank pressures
and setpoints and saves/prints total wall time including setup, solve, and plots.

### 7.4 Reproduced early nonlinear stall

The stall was independently reproduced with regen disabled at the attempted
step from 0.015 to 0.020 simulated seconds. FullFlow's default SciPy nonlinear
solve could spend thousands of evaluations on a difficult attempt before
returning control to the outer timestep-retry loop.

`Hopper_Solver.py` now subclasses the installed FullFlow `TransientStepSolve`
and supplies SciPy `max_nfev` for each least-squares attempt. Exhausting the
budget still goes through FullFlow's residual check and existing smaller-step
retry. A failed attempt is not accepted because its budget expired.

The tested 80-evaluation cap crossed the initial difficult step by retrying at
2.5 ms. The final default is **20**, exposed as
`simulation.nonlinear_max_evaluations`; older configs get that materialized default.
The maximum step remains 5 ms for these inputs and accepted residual tolerance
remains `1e-6`.

This caps the least-squares main evaluation count; numerical Jacobian probes
are additional work. It is not a strict iteration count, total-residual-call
limit, whole-solver wall-clock timeout, or five-minute deadline. It also does
not cap every internal Newton/Broyden path in FullFlow. The adapter uses private
FullFlow interfaces, so dependency upgrades need focused compatibility checks.
No installed package source was patched.

### 7.5 Performance and numerical comparison

Two full cold-start driver runs of the corrected **10-second ramp / 11-second
sequence**, including all plots, measured:

- Cap 80: 439.6 s total; 2,230 accepted steps; 31 retries;
  maximum accepted residual below `1e-6`; gas mass drift below 9 mg.
- Cap 20: 298.8 s total; 2,264 accepted steps; 123 retries;
  maximum accepted residual below `1e-6`; gas mass drift below 9 mg.
- Demonstrated maximum: 2055.91 versus 2062.31 N, a 0.31% change.
- Final LOX/fuel inventory differences: 5.42 g / 3.03 g.
- Maximum pointwise tank-pressure difference: about 20.7 psi, associated with
  changed bang-bang pulse timing. A fixed residual tolerance does not prove
  identical transient trajectories.

The five-minute target was met narrowly in the tuning benchmark, but a later
committed 10-second run took **310.96 s**. Therefore **a reliable sub-five-minute
guarantee has not been established**. Input differences, retries, dependency
versions, CPU load, and property work affect runtime. Shortening a ramp is not
evidence that the original full profile became faster.

The performance benchmark directories/scripts are local ignored development
artifacts under `validation_outputs/`, not portable fixtures delivered by Git.
The committed run directories provide additional, independently saved evidence.

### 7.6 Interpreter/launcher problem

The default terminal `python` resolved to an MSYS Python without matplotlib.
The old `venv` pointed to a missing Python 3.11 installation. Separately, the
user's Git Bash invocation contained a malformed path starting `:/Users/...`,
which failed before Python started.

A working `.venv-validation` environment was created locally and used for tests.
It is ignored and is not part of the portable repository. `run_hopper.sh` changes
to its own directory and selects a working Windows environment from `.venv`,
`.venv-validation`, then `venv`, checking the required imports before execution.
It avoids editor-generated drive-prefix paths and forwards CLI arguments.
It is a Git Bash/Windows launcher, not a generic Unix environment bootstrap.

## 8. Final checked-in inputs

These are the values in `996f8509`, not timeless hardware specifications and
not necessarily the values used in every earlier benchmark. Consult JSON and
each run's snapshots before comparing results.

### 8.1 Simulation, sequence, and engine

- Selected profile: linear ramp, 100 to 600 lbf over **2.5 s**, **40 points**.
- Requested maximum step: 0.05 s; save interval: 0.05 s; resolved step cap: 0.005 s.
- Residual tolerance: `1e-6`; Jacobian: `2-point`; nonlinear cap: 20.
- Controller period: 0.01 s; two steps per period; ten steps per valve full stroke.
- Initial valve-demand start: 0.3 s; demand ramp: 0.2 s; tail-off: 0.3 s.
  With this ramp the total nominal sequence is 3.5 simulated seconds.
- Target O/F: 1.25. Throat area: 1.0420050840101682 in²; exit area:
  3.651185814371629 in²; area ratio approximately 3.504; throat diameter about
  1.1518 in.
- Ambient pressure: 14.8 psia; nominal chamber pressure: 300 psia.
- Map: 10–610 psia, 31 points; MR 0.5–4.0, 36 points.
- C-star efficiency: 0.80 to 0.85 from 100 to 300 psia, held at edges.
- Cf efficiency: 0.94 to 0.97 from 100 to 250 psia, held at edges.
- Sustained-thrust hold: 0.5 s; tracking shortfall criterion: 2%.
- Property-temperature table: 101 points; flow floor `1e-4` kg/s;
  nozzle small-drop regularization 0.5 psi.
- All plots and diagram enabled; plot CSV output disabled in the final config.

### 8.2 Sources and tanks

- Fluids: LOX, RP-1, and nitrogen pressurant.
- Two 9 L bottles, each at 4,503.2 psia and 293 K; each temperature limit 220 K;
  outlet CdA 100 mm². Common-header volume 0.05 L.
- LOX tank: 20 L; initial mass 19.43 kg; liquid 90 K; initial ullage 150 K;
  initial/set pressure 500 psia; deadbands -5/+5 psi; surface radius about
  7.874 in; HTC 500 W/(m² K).
- Fuel tank: 16 L; initial mass 9.7 kg; liquid/ullage 293 K; initial/set pressure
  500 psia; deadbands -5/+5 psi; surface radius about 3.937 in; HTC zero.
- Both initial controller commands are one; physical gas solenoid travel starts
  from closed in the builder.

### 8.3 Ordered flow paths

- Oxidizer pressurant: 48 in tube, OD 0.25 in/wall 0.035 in, inertia disabled;
  OBB CdA 18.68 mm² with 0.05/0.05 s travel; OR-1 CdA **6.04 mm²**;
  PROCK CdA 12 mm² with ideal timing.
- Fuel pressurant: matching tube dimensions; FBB CdA 18.68 mm² with
  0.05/0.05 s travel; OR-2 CdA **6.02 mm²**; PRFCK CdA 12 mm² with ideal timing.
- Liquid tubing: OD 0.375 in/wall 0.035 in, giving ID 0.305 in.
- LOX feed: 48 in tube with zero height change; OTV max CdA **50 mm²** and
  0.2/0.2 s travel; 10 in tube; OINJ CdA 15.3 mm² and 0.3937 in passage.
- Fuel feed: 12 in tube; FTV max CdA **50 mm²** and 0.2/0.2 s travel;
  12 in tube; regen CdA 25 mm² and 11.811 in passage; FINJ CdA 10.6 mm² and
  0.3937 in passage.
- Regen curve: `[150 psia, 51 K rise]`, `[300 psia, 44 K rise]`.
- Tube roughness: about 0.00059055 in.

The early migration deliberately kept gas-line OD/wall unspecified because
only ID had been supplied then. The final committed config now supplies
quarter-inch OD and 0.035-inch wall. Older README/DESIGN statements about
unspecified gas dimensions or zero default regen describe earlier stages;
the actual config/run snapshot takes precedence.

## 9. File-by-file implementation guide

The manifest supplies exact per-commit file statuses. This section explains
why the functional files matter.

### 9.1 Driver and configuration

- [Hopper.py](Hopper.py): CLI config selection/results-only switch, environment
  setup, initial-pressure output, archive-before-import order, map/target/solve
  orchestration, diagram/reports/exports, and saved total runtime.
- [hopper_config.json](hopper_config.json): complete public vehicle/run inputs
  for the serialized built-in model, including user parameter edits.
- [HOPPER_DATA.py](HOPPER_DATA.py): reduced from editable Python inputs to a
  compatibility adapter, preserving raw source text and exporting SI values.
- [Hopper_Config.py](Hopper_Config.py): JSON reading, duplicate/unknown-key checks,
  validation, object construction, efficiency modes, profile selection,
  derived values, and integration-step calculation.
- [Hopper_Config_Fields.py](Hopper_Config_Fields.py): public field names,
  SI conversion factors, and legacy-symbol mappings.
- [Hopper_Run_Config.py](Hopper_Run_Config.py): unique archive creation, input
  asset copying, portable paths, hashes, and derived input snapshots.
- [Hopper_Inputs.py](Hopper_Inputs.py): shared ramp/profile checks, cached CSV
  reading/interpolation, efficiency/rise validation, and asymmetric slew math.

### 9.2 Physics and numerical behavior

- [Hopper_Elements.py](Hopper_Elements.py): tank/source data objects, OD/wall/ID
  validation, line geometry, closed-safe inertial restrictions, valve actuator
  wiring, demand-versus-physical throttle separation, Regen, COPVSource,
  and gas/check-valve handling.
- [Hopper_Control.py](Hopper_Control.py): accepted-state helper, dictionary
  FunctionLookup correction, retry-safe ValveTravel, sampled hysteresis and
  controller breakpoints.
- [Hopper_Prop_System.py](Hopper_Prop_System.py): individual bottle/header
  construction, tank and ordered path builders, per-component gas openings,
  heated liquid contexts/property table, additional tracks, step cap,
  bounded solver use, startup progress, and post-solve metadata writing.
- [Hopper_Components.py](Hopper_Components.py): shared engine efficiency lookup
  applied consistently to targeter and transient; existing nozzle/map block.
- [Hopper_Thrust_Targeter.py](Hopper_Thrust_Targeter.py): profile validation/sorting
  before steady Pc and flow-demand solves.
- [Hopper_Solver.py](Hopper_Solver.py): FullFlow adapter that adds `max_nfev`
  while retaining residual acceptance and retry behavior.
- [Hopper_Chamber_Map.py](Hopper_Chamber_Map.py) and
  [Hopper_Chamber_Map.h5](Hopper_Chamber_Map.h5): existing equilibrium gas-map
  generator/data. They were moved back from the remote `Backend/` layout in
  the replacement commit; the ideal-map concept was not replaced by an
  empirical efficiency map.

### 9.3 Output, documentation, and support

- [Hopper_Export.py](Hopper_Export.py): existing target/actual-performance CSV
  exports, extended with optional thermal tracks.
- [Hopper_Report.py](Hopper_Report.py): HDF5 metadata/config, sustained capacity,
  inventory and temperature-limit summaries, validity flags, and report files.
- [Hopper_Plots.py](Hopper_Plots.py): station snapshots, regen/COPV plots,
  heated-property histories, fractional gas-drop allocation, hardware/path
  checks, string CSV handling, and removal of Hdot/quality panels.
- [Hopper_Diagram.py](Hopper_Diagram.py): configuration-derived PNG/SVG schematic,
  compact branch routing, symbols, imperial geometry, and derived throat diameter.
- [run_hopper.sh](run_hopper.sh): checked Windows-environment selection for Git Bash.
- [requirements.txt](requirements.txt): existing dependencies; FullFlow pinned
  to 3.0.0, with several other packages unpinned or lower-bounded.
- [.gitignore](.gitignore): environment/cache/scratch/raw-history exclusions.
- [README.md](README.md): run/config/plot/reference workflow and limits.
- [DESIGN.md](DESIGN.md): progress checklist, decisions, and validation history.
- [inputs/cstar_efficiency.example.csv](inputs/cstar_efficiency.example.csv) and
  [inputs/cf_efficiency.example.csv](inputs/cf_efficiency.example.csv): format
  examples, not measured engine calibration.
- [Hopper_Results/](Hopper_Results/): tracked sample outputs and three archived
  run input/report sets; exact files are in the manifest.
- [comparison/](comparison/): historical comparison scripts and figures. Their
  previous remote `Old_sims_comparison`/`Outputs` layout was reverted by the
  folder replacement. Historical comparisons were not regenerated as a new
  validation baseline for all these physics changes.
- [thrust_profile_10m_tethered.csv](thrust_profile_10m_tethered.csv) and
  [thrust_profile_52m_free.csv](thrust_profile_52m_free.csv): user trajectory data
  included in the second commit; not automatically consumed by the loader.

### 9.4 Validation code

- [tests/test_features.py](tests/test_features.py): 12 tests for live function
  lookup, tube geometry, CSV interpolation/errors, ramp/rise behavior,
  actuator travel/retry safety, controller hysteresis, and report behavior.
- [tests/test_config.py](tests/test_config.py): seven tests covering config
  conversions/defaults, invalid/duplicate fields, geometry/timing rejection,
  efficiency modes/relative paths, independent archives, and bounded-solver
  defaults/validation for old configs.
- [tests/test_solver.py](tests/test_solver.py): two tests checking invalid
  budgets and that the adapter supplies the work budget without relaxing
  acceptance-related settings. These are interface tests, not a complete
  nonlinear-solver correctness proof.
- [tests/run_feature_validation.py](tests/run_feature_validation.py): isolated
  in-process 500-to-700 N short coupled run, nonzero regen and a deliberately
  violated COPV limit, reports, integration check, all plots/CSV output.
- [tests/check_integration.py](tests/check_integration.py): saved-run checks
  for heating, downstream temperature, valve rates, inventory conservation,
  bottle behavior, reports, property interpolation, and performance CSV export.
  Bottle names now come from saved metadata rather than a hardcoded count of three.

## 10. Validation and its limits

### 10.1 Unit checks

The suite grew from 12 to 18 and finally **21 passing tests**. It was rerun
against the exact replacement tree before publication and again while preparing
this handoff. No simulator behavior was changed for this documentation task.

### 10.2 Coupled and output checks performed during development

- A complete short 1.5 s startup/burn/shutdown feature run at 5 ms integration
  and output intervals exercised a 500-to-700 N request, nonzero regen, and a
  deliberate bottle temperature-limit violation.
- All 19 named plots rendered with plot CSVs; station, regen, and schematic
  outputs were visually inspected.
- Tank mass consumption was checked against integrated tank outflow. Total
  connected gas inventory was used for pressurant conservation checks.
- At discrete valve events, saved connector rates are endpoint rates, not
  interval-averaged transfers. A per-bottle endpoint-flow integral did not
  perfectly match inventory; that distinction is documented rather than
  presenting it as a whole-system mass loss.
- Checked heated-property interpolation errors were below 0.003% at sampled
  points for the tested range. The checker accepts a broader `rtol=2e-3`;
  neither statement proves accuracy over arbitrary user ranges.
- A JSON-driven short run verified copied assets, asset hashes, embedded
  config, all plots/CSVs, and results-only re-export without rerunning the
  transient.
- The cold 10-second benchmark and 80-versus-20 budget comparison are described
  in section 7. They did not establish temporal convergence of pressure pulses.

The short feature runner depends on recognizable component names/timings and
the active hardware configuration. Its input overrides are in-process, so it
is not a completely frozen hardware fixture or a faithful example of a driver
input archive. Use an explicit saved case when validating a changed design.

### 10.3 Committed run evidence

The following data were read directly from committed `run_config.json`,
`run_timing.json`, and summaries while preparing this handoff:

- [20261105T095434_244612Z](Hopper_Results/20261105T095434_244612Z/): 10.0 s,
  101-point ramp; **310.96 s** total; sustained demonstrated maximum about
  **2062.3 N / 463.6 lbf**; final LOX/fuel **14.2985 / 6.0136 kg**; bottle
  minimum **223.25 K**.
- [20261105T100133_895127Z](Hopper_Results/20261105T100133_895127Z/): 4.0 s,
  40-point ramp; **281.92 s** total; demonstrated maximum about
  **2032.0 N / 456.8 lbf**; final LOX/fuel **17.2753 / 8.1430 kg**; bottle
  minimum **269.59 K**.
- [20261105T100916_169935Z](Hopper_Results/20261105T100916_169935Z/): 2.5 s,
  40-point ramp; **212.00 s** total; demonstrated maximum about
  **2050.3 N / 460.9 lbf**; final LOX/fuel **18.0204 / 8.6750 kg**; bottle
  minimum **277.88 K**.

All three summaries report sustained valve saturation/shortfall evidence and
no bottle limit violation at the saved 220 K threshold. They have different
input histories and should not be treated as a controlled duration-only study.
Reported numbers describe simulated cases; they are not a hardware performance
certificate or a universal maximum-thrust rating.

## 11. Running, replaying, and extending the code

### 11.1 Environment and ordinary run

Run from `sims/hopper_fullflow` with a working Windows CPython environment and
the packages in `requirements.txt`. The local validated environment used Python
3.12 with FullFlow 3.0.0, FullPlot 2.0.0, and ThermoProp 2.1.2. The old venv path
or the shell's default `python` may select a different/broken interpreter.

For the existing local environment:

```powershell
.\.venv-validation\Scripts\python.exe -m unittest discover -s tests -v
.\.venv-validation\Scripts\python.exe Hopper.py
.\.venv-validation\Scripts\python.exe Hopper.py --config path/to/case.json
```

From Git Bash:

```bash
bash run_hopper.sh
bash run_hopper.sh --config path/to/case.json
```

On a fresh checkout, create a Windows Python environment and install
`requirements.txt`; the ignored local validation environment is not provided
by Git. When using a different interpreter, verify package imports first.
Do not reinstall into an unrelated MSYS Python merely because it is first on PATH.

### 11.2 Replay results or rerun a case

```powershell
.\.venv-validation\Scripts\python.exe Hopper.py --config Hopper_Results/<run>/run_config.json --results-only
```

This requires the raw HDF5 history locally. It does not rerun the transient,
but the driver still performs its map check and thrust-target calculations.
Use an existing map and disable rebuild for a pure re-export workflow; the
results-only switch does not independently suppress every upstream operation.
Re-export can replace output files in that archived directory, so preserve a
copy if exact original figures/reports must remain immutable.

To rerun the archived numerical case, use the same config without
`--results-only`; the normal new-run archive workflow creates a new timestamped
subdirectory under its configured results root. For important work, preserve
the matching code/dependency versions too.

For standalone plots, select the run's config before importing plotting modules:

```powershell
$env:HOPPER_CONFIG = 'Hopper_Results/<run>/run_config.json'
.\.venv-validation\Scripts\python.exe Hopper_Plots.py list
.\.venv-validation\Scripts\python.exe Hopper_Plots.py plot ox_pressure_stations fuel_pressure_stations --times 1 2 --csv
Remove-Item Env:HOPPER_CONFIG
```

Times must be relevant to the selected run. `--profile-time` changes the plot
clock to profile start; default history time includes startup sequencing.

### 11.3 Extend inputs and physics deliberately

- Add/reorder supported component specs in JSON to change hardware topology;
  inspect the regenerated diagram and station labels before trusting results.
- For a new component class, implement the builder in `Hopper_Elements.py`,
  add explicit serialized fields/unit conversions and type support in the config
  loader, and extend tracking/report tests where behavior changes.
- Keep efficiency logic shared between targeter and transient.
- Do not introduce calculated dictionary lookups without preserving the
  `FunctionLookup` reevaluation behavior.
- Keep controllers/actuators based on accepted states so rejected attempts do
  not advance control state.
- Regenerate maps when reactants/thermodynamic grid change, and validate the
  intended property/model range.
- For changes to startup, gas storage, solver budget, or sample period, check
  time histories and conservation as well as completion/runtime.

## 12. Known limitations and follow-up priorities

### 12.1 Highest priority: numerical and physical confidence

1. **Temporal convergence of tank-pressure pulses.** Repeat a fixed, archived
   hardware case with progressively smaller integration and saved-output steps.
   Compare pulse envelopes/timing, extrema, inventories, delivered thrust, and
   operating-limit crossings. The 20.7 psi comparison difference remains an
   unresolved accuracy question, not a solved convergence result.
2. **Reliable runtime margin.** Profile accepted and rejected nonlinear work,
   investigate the stiff small header/large bottle connections, and consider
   caching property-table construction or better solver initialization. The
   current measured five-minute result has little margin; a committed longer
   result exceeds it. Avoid coarsening time or relaxing tolerance solely to
   claim a speed target.
3. **Hardware calibration.** Supply measured bottle/header plumbing, valve
   travel and CdAs, tank thermal parameters, regen curves, efficiencies, and
   actual line geometry. Placeholder values are not validated hardware data.
4. **Regen/combustion coupling.** Decide whether to add wall thermal states,
   phase validity checks, and reactant enthalpy feedback into the combustion
   model. Present heating is prescribed and tank liquid stays isothermal.

### 12.2 Reproducibility and software follow-up

- Freeze dependency versions for repeatable engineering runs. Only FullFlow
  is exactly pinned; private solver APIs and callback behavior can change.
- Add self-contained integration fixtures independent of the editable vehicle
  config, with explicit assertions for the difficult early timestep.
- Extend whole-system gas conservation checks for dynamic gas-line configurations;
  the present check sums bottles, header, and tank ullages and assumes the
  default lumped gas legs. Additional stored gas junctions need inclusion.
- Consider input JSON Schema/editor help and more complete physical/type
  validation. Add direct CSV-thrust loading if desired; the trajectory files
  included here are currently data only.
- Add automatic map/config compatibility checks rather than relying only on
  file existence and user-directed rebuilds.
- Audit header initialization for heterogeneous bottles: the network initializes
  header pressure from minimum bottle pressure and temperature from the mean,
  while `resolved_inputs.json` currently reports compatibility constants based
  on the first bottle. Identical checked bottles hide this provenance mismatch.
  Correct the saved derived metadata before relying on mixed initial conditions.
- Expand archived derived values and record code/environment identifiers if
  exact long-term replay is required. Keep raw HDF5 histories separately.
- Refactor the top-level driver into a library-safe `main()` if embedding,
  repeated cases in one process, or automated parameter sweeps are needed.
- Port remote-only Monte Carlo/isolated-analysis features only after deciding
  their desired relationship to this JSON configuration and model. They were
  intentionally overwritten, not automatically retained.

### 12.3 Physical/model boundaries that still apply

The chamber is quasi-steady and zero-volume; ignition is immediate when flow
arrives. Sensors, reliefs, purge circuits, vents, drains, ground interfaces,
line priming, and other omitted hardware are not made present merely by drawing
the existing component graph. Feed junctions have no trapped liquid inventory.
Gas intermediate station pressures in lumped groups are approximations. COPV
temperature limits are flags. Tank depletion is not a two-phase simulation.
The five-minute goal, an all-hardware capacity guarantee, full pressure-pulse
convergence, and experimental calibration remain separate acceptance questions.

## Appendix: audit the changes from Git

From the repository root, these commands reproduce the per-commit source view:

```bash
git show --stat 31d3aad8 -- sims/hopper_fullflow
git show --stat 9b85078b -- sims/hopper_fullflow
git show --stat 996f8509 -- sims/hopper_fullflow
git diff 31d3aad8^ 31d3aad8 -- sims/hopper_fullflow
git diff 9b85078b^ 9b85078b -- sims/hopper_fullflow
git diff 996f8509^ 996f8509 -- sims/hopper_fullflow
```

The companion [file manifest](docs/handoff_commit_files.txt) captures all
three complete path/status lists. It includes unchanged-content renames as Git
detects them, binary figure/map changes, user-provided trajectory data, sample
results, and the remote-only files removed by the final folder replacement.
This documentation addition is subsequent to those three implementation commits.

# Test-Stand Hotfire Sim: Inputs

Sim: `isolated_sims/Hopper_Hotfire.py` · Inputs: `isolated_sims/hotfire_config.json` · Fluids, chamber map, ambient pressure and valve sequence timing: `hopper_config.json`
Pressures are absolute (psia) unless marked psig (gauge) or psi (difference). **GUESS** = nobody has confirmed it. **COPIED** = copied from `hopper_config.json`.

```
python isolated_sims/Hopper_Hotfire.py                 # run with isolated_sims/hotfire_config.json
python isolated_sims/Hopper_Hotfire.py --check-inputs  # build the hardware and print the timeline only
python isolated_sims/Hopper_Hotfire.py --config my_case.json --out some/folder
```

## Layout

Each side runs tank → test-stand CdA → throttle valve → line → injector → chamber. On the fuel side the regen jacket sits between the line and the injector.

```
LOX tank  -> LOX Stand  -> OTV -> LOX Line  ----------------> OINJ -> Main Chamber -> nozzle
Fuel tank -> Fuel Stand -> FTV -> Fuel Line -> Regen Jacket -> FINJ -> Main Chamber
```

The torch igniter and the purges are not modelled.

## Reservoirs, stand and lines (per side)

| Input | LOX | Fuel | Notes |
|---|---|---|---|
| Reservoir set pressure | 800 psia | 800 psia | Held constant: no blowdown, pressurant or ullage |
| Liquid temperature | 90 K | 293 K | COPIED |
| Test-stand CdA | 33 mm² | 28 mm² | All stand plumbing upstream of the throttle valve as one restriction |
| Throttle valve full-open CdA | 24 mm² | 24 mm² | COPIED |
| Throttle valve actuation | 0.18 s | 0.18 s | COPIED |
| Stand-to-injector line | 72 in, 3/8 in OD × 0.028 in wall | 72 in, 3/8 in OD × 0.028 in wall | Roughness 0.00059 in |
| Regen | n/a | 30 mm², +82 K at 150 and 300 psia Pc | COPIED |
| Injector CdA | 16.6 mm² (OINJ) | 15.0 mm² (FINJ) | COPIED |

## Engine

These fields are the same as the `engine` section of `hopper_config.json`.

| Input | Value | Notes |
|---|---|---|
| Mixture ratio | 1.25 | Regen TCA CDR |
| Throat area | 1.0568 in² | Regen TCA CDR |
| Exit area | 3.0946 in² | CDR exit radius 0.9925 in; expansion ratio 2.928 |
| Nominal Pc | 300 psia | Regen TCA CDR |
| c* efficiency | 0.80 at 100 psia → 0.85 at 300 psia | Linear in Pc, flat outside, as in `hopper_config.json`. 0.85 is the CDR value |
| Cf efficiency | 0.94 at 100 psia → 0.97 at 250 psia | Linear in Pc, flat outside. COPIED (the CDR gives none) |

## Thrust profile

| Input | Value | Notes |
|---|---|---|
| CSV | `thrust_curves/hotfire_ramp.csv` | Columns `time [s]`, `thrust [lbf]`; path relative to the hotfire config |
| Profile | 0 → 550 lbf over 5 s (a point every 0.1 s), hold 550 lbf for 1 s | User, 2026-11-19 |
| Interpolation | linear | `linear` or `step` |
| 0 lbf points | Both throttle valves shut | Zero flow targets; the targeter solves every other point |
| Profile time 0 | 0.7 s on the sim clock | Main ignition (`sequence.main_ignition_s`); the last point is the shutdown request |

## Simulation

| Input | Value | Notes |
|---|---|---|
| Maximum time step | 0.015 s | The solver also lands on every valve and profile breakpoint |
| Save interval | every step | |
| Output folder | `Outputs/Hotfire/<UTC time>/` | |

## Outputs

| File | Contents |
|---|---|
| `hotfire_history.csv` | Every time step: thrust, target, Pc, MR, Isp, c*/Cf efficiency, exit pressure; for each side the tank pressure and mass flow; each element's inlet pressure and dP; throttle valve target flow, CdA and % open; injector stiffness; regen dT |
| `hotfire_summary.json` | Peak thrust and the state there; when each throttle valve first hits full open; minimum stiffness and TV dP; maximum Pc, TV dP, % open and injector dP; settled values for each hold; propellant used |
| `hotfire_holds.csv` | Settled values (second half of each constant-thrust hold) |
| `thrust_targets.csv` | Targeter Pc and flows for each profile point |
| `hotfire_overview.png` | Time histories: thrust, Pc, mass flows, MR, TV % open, TV CdA, TV dP, injector dP, stiffness |
| `hotfire_vs_thrust.png` | Pc, TV % open, TV dP, stiffness, mass flows and Isp against thrust |
| `hotfire_pc_vs_thrust.png` | Chamber pressure against thrust |
| `hotfire_element_dp.png` | dP of every element on each side over time |
| `hotfire_pressure_ladder.png` | Pressure at each station at peak thrust |
| `hopper_config_hotfire_engine.json` | The `hopper_config.json` the run used, with the hotfire engine in it |
| `hotfire_config.json`, profile CSV | Copies of the inputs; the copied config points at the copied CSV |

## Modelling assumptions

- Run tanks are constant-pressure reservoirs.
- The test-stand CdA is one fixed incompressible restriction. It has a 1 cm passage length for flow inertia only, which adds no steady dP.
- OTV/FTV deliver the thrust-targeter flows up to their full-open CdA. Past that the valve stays fully open and the thrust falls short of the target.
- Propellant ignites as soon as it reaches the chamber. The chamber has no volume.
- Liquids are incompressible at the reservoir pressure and temperature. The regen heats the fuel by the configured dT.

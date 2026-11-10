# Test-Stand Hotfire Sim: Inputs

Sim: `isolated_sims/Hopper_Hotfire.py` · Inputs: `isolated_sims/hotfire_config.json` · Engine data: `hopper_config.json`
Pressures are absolute (psia) unless marked psig (gauge) or psi (difference). **GUESS** = nobody has confirmed it. **COPIED** = copied from the vehicle config, not measured on the stand.

## Simulation

| Input | Value | Notes |
|---|---|---|
| Maximum time step | 0.02 s | Solver also lands on every valve and thrust-ramp breakpoint |
| Save interval | every step | |
| Tail-off after last event | 0.3 s | |
| Output folder | `Outputs/Hotfire/<UTC time>/` | |

## Run tanks (reservoirs)

| Tank | Propellant | Pressure | Temperature | Notes |
|---|---|---|---|---|
| LOX Tank | LOX | 1000 psia | 90 K | Constant pressure: no blowdown, pressurant or ullage |
| Kerosene Tank | RP-1 | 1000 psia | 293 K | Constant pressure |

## Thrust profile

| Input | Value | Notes |
|---|---|---|
| Start thrust | 100 lbf | |
| End thrust | 500 lbf | |
| Ramp duration | 2.5 s | |
| Ramp points | 26 | One every 0.1 s |
| Hold at end thrust | 0 s | Set > 0 to add a steady hold |
| Interpolation | linear | `linear` or `step` |
| Profile time 0 | 0.5 s on the sim clock | When OTV/FTV finish opening (from `hopper_config.json`) |

## Engine (from `hopper_config.json`)

| Input | Value | Notes |
|---|---|---|
| Target mixture ratio | 1.25 | |
| Throat area | 1.042 in² | |
| Exit area | 3.651 in² | Expansion ratio 3.50 |
| c* efficiency | 0.80 at 100 psia → 0.85 at 300 psia | Linear, flat outside |
| Cf efficiency | 0.94 at 100 psia → 0.97 at 250 psia | Linear, flat outside |
| Ambient pressure | 14.8 psia | |
| OTV/FTV start opening | 0.3 s | |
| OTV/FTV opening (and closing) time | 0.2 s | |

## LOX feed

| Component | Type | Value | Notes |
|---|---|---|---|
| LOX Tank to OTV Line | Line | 3/8 in OD × 0.010 in wall (0.355 in ID), 48 in long | COPIED |
| OTV | Throttle valve | 50 mm² full open | Opens to whatever CdA meets the flow target. COPIED |
| OTV to LOX Injector Line | Line | 3/8 in OD × 0.010 in wall, 10 in long | COPIED |
| OINJ | Injector | 16.6 mm² | COPIED |

## Fuel feed

| Component | Type | Value | Notes |
|---|---|---|---|
| Fuel Tank to FTV Line | Line | 3/8 in OD × 0.010 in wall, 12 in long | COPIED |
| FTV | Throttle valve | 50 mm² full open | COPIED |
| FTV to Regen Line | Line | 3/8 in OD × 0.010 in wall, 12 in long | FRMPT reads this manifold. COPIED |
| Regen Jacket | Regen | 30 mm²; +82 K at 150 and 300 psia Pc | COPIED |
| FINJ | Injector | 15.0 mm² | COPIED |

All lines: roughness 0.00059 in, no height change.

## Torch igniter

| Input | Value | Notes |
|---|---|---|
| OIMV / FIMV open | 0.2 s | |
| OIMV / FIMV close | 2.0 s | |
| OIMV / FIMV stroke | 0.08 s | |
| Torch throat area | 0.006504 in² | About 0.091 in diameter |
| Torch c* efficiency | 1.0 | GUESS (ideal) |
| Zero-flow smoothing band | 0.5 psi | Numerical only |

**Torch oxidizer leg** (tees off the LOX line above OTV)

| Component | Type | Value | Notes |
|---|---|---|---|
| OIMV | Valve | 0.68 mm² | |
| OICK | Check valve | 5.0 mm² | |
| OIORF | Orifice | 8.0 mm² | Provisional |
| OTHX | Tube | 1/4 in OD × 0.065 in wall (0.120 in ID), 60 in long | Treated as a plain tube: no heating or vaporization |
| Torch Ox Injector | Injector | 0.4398 mm² | |

**Torch fuel leg** (tees off the kerosene line above FTV)

| Component | Type | Value | Notes |
|---|---|---|---|
| FIMV | Valve | 2.2 mm² | |
| FICK | Check valve | 5.0 mm² | |
| FIORF | Orifice | 8.0 mm² | Provisional |
| Torch Fuel Injector | Injector | 0.8051 mm² | |

## Purge system

| Input | Value | Notes |
|---|---|---|
| GN2 supply pressure | 2000 psia | GUESS, held constant |
| GN2 supply temperature | 293 K | |
| Purge Reg set pressure | 750 psig | Ideal regulator (no droop). From PR-REG |
| Purge Reg CdA | 5.0 mm² | GUESS. From PR-REG |
| FIPV opens | 2.0 s | When the torch valves close |
| FIPV closes | 1.0 s after OTV/FTV are shut | GUESS |
| OPV / FPV open | 0 s after OTV/FTV are shut | GUESS |
| OPV / FPV open for | 1.0 s | GUESS |
| Purge valve stroke | 0.01 s | |

**Purge paths** (all start at the Purge Reg outlet)

| Path | Components | Ends at | Notes |
|---|---|---|---|
| Purge Relief | Purge RV: 2.0 mm², cracks at 825 psig | Ambient | From ORV-1. GUESS |
| LOX Purge | Line 1/4 in × 0.035 in wall (0.180 in ID), 39.4 in long → OPV 1.86 mm² → OPCK 8.5 mm² | LOX line below OTV | Line length GUESS; OPCK GUESS |
| Fuel Purge | Line 1/4 in × 0.035 in wall, 39.4 in long → FPV 1.86 mm² → FPCK 8.5 mm² | Fuel line below FTV | Line length GUESS; FPCK GUESS |
| Torch Purge | Line 1/4 in × 0.035 in wall, 39.4 in long → FIPV 2.2 mm² → FIPCK 5.0 mm² → FIPORF 50 mm² → torch fuel injector 0.8051 mm² | Torch chamber | FIPV is called TFPV in `hopper_config.json`. Keep the injector CdA equal to the torch fuel leg's |

## Modelling assumptions

- Run tanks are constant-pressure reservoirs.
- OTV/FTV deliver the thrust-targeter flows. The targets ignore the torch flow.
- Propellant ignites as soon as it reaches a chamber. Neither chamber has volume.
- The torch burns on the same LOX/RP-1 chamber data as the main engine. Its throat flow is ideal and can run either way.
- Liquids are incompressible at the reservoir pressure and temperature.
- Purge gas only sees its back pressure. It is not added to the propellant lines or chambers. Check valves block reverse flow.
- The torch purge is solved as if FIMV is shut.

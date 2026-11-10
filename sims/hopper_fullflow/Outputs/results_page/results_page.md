# Free-Flight Profile: 52m_free (full vehicle)

[image: flight_overview.png]

[image: startup_sequence.png]

[image: shutdown_sequence.png]

[file: free flight simulation.csv]

| Graph | Quantity | Max |
| --- | --- | --- |
| Mass flow | LOX, main injector | 0.5636 kg/s |
|  | Fuel, main injector | 0.4509 kg/s |
|  | LOX, torch | 8.72 g/s |
|  | Fuel, torch | 8.64 g/s |
| Throttle valves | OTV open | 25.2 % |
|  | FTV open | 26.8 % |
| Thrust | Thrust | 379.7 lbf |
|  | Target thrust | 379.7 lbf |
| Chamber pressure | Main chamber Pc | 273.7 psia |
|  | Torch chamber Pc | 552.6 psia |
|  | Torch over main chamber | 537.7 psi |
| Injector dP | LOX (OINJ) | 72.7 psi |
|  | Fuel (FINJ) | 87.0 psi |
| Mixture ratio | Main chamber MR | 4 |
|  | Torch MR | 4 |
| Tank pressure | LOX tank | 566.0 psia |
|  | Fuel tank | 569.6 psia |
| Propellant load | LOX liquid | 16.0 kg |
|  | Fuel liquid | 13.0 kg |
| Ullage temperature | LOX ullage | 206.7 K |
|  | Fuel ullage | 300.4 K |
| COPV | COPV pressure | 4503.2 psia |
|  | COPV gas temperature | 293.0 K |
|  | COPV gas mass | 5.58 kg |
| Pressurant flow | LOX pressurant (OBB) | 314.1 g/s |
|  | Fuel pressurant (FBB) | 198.7 g/s |
| Purge regulator | PR-REG outlet | 819.1 psia |
|  | PR-REG opening | 100.0 % |
|  | PR-REG flow | 24.3 g/s |
| Purge flow | LOX purge (OPV) | 24.0 g/s |
|  | Fuel purge (FPV) | 14.3 g/s |
|  | Torch purge (TFPV) | 12.8 g/s |
| Regen | Regen outlet temperature | 375.0 K |
|  | Regen jacket dP | 20.2 psi |
| Relief flow | ORV (LOX tank) | 9.45 g/s |
|  | FRV (fuel tank) | 11.8 g/s |
|  | PRRV (COPV) | 0 g/s |
|  | ORV-1 (purge reg) | 0 g/s |

## Run summary

| Metric | Value |
| --- | --- |
| Burn time [s] | 29.85 |
| Total impulse [N s] | 38890.0 |
| LOX residual [kg] | 2.279 |
| LOX residual [%] | 14.25 |
| Fuel residual [kg] | 1.96 |
| Fuel residual [%] | 15.08 |
| LOX used [kg] | 13.719 |
| Fuel used [kg] | 11.04 |
| COPV pressure at shutdown [psia] | 1346.6 |
| COPV end pressure [psia] | 1311.3 |
| COPV gas used [kg] | 2.54 |
| COPV gas left [kg] | 3.045 |
| LOX tank min in profile [psia] | 489.9 |
| Fuel tank min in profile [psia] | 488.6 |
| LOX tank max [psia] | 566.0 |
| Fuel tank max [psia] | 569.6 |
| Min LOX injector stiffness [% Pc] | 19.1 |
| Min fuel injector stiffness [% Pc] | 22.8 |
| Max /thrust error/ at profile points [%] | 1.25 |
| RMS thrust error at profile points [%] | 0.08 |
| Profile points judged | 295 of 296 |
| Max OTV % of max CdA | 25.2 |
| Max FTV % of max CdA | 26.3 |
| Relief mass vented [kg] | 0.0016 |
| All thrust points reachable | yes |
| Thrust within tolerance | yes |
| Max OTV CdA needed at current set pressure [mm2] | 12.41 |
| Max FTV CdA needed at current set pressure [mm2] | 13.08 |
| COPV end-of-life pressure [psia] | 1380.1 |
| COPV end-of-life temperature [K] | 209.6 |
| Regulator outlet [psia] | 764.8 |
| Regulator outlet temperature [K] | 195.9 |
| Regulator capacity at EOL [g/s] | 24.03 |
| LOX purge path CdA [mm2] | 1.75 |
| LOX max propellant volumetric flow [L/s] | 0.49 |
| LOX purge flow at that moment [L/s] | 0.503 |
| LOX purge / propellant at peak [%] | 102.5 |
| LOX feed pressure at peak [psia] | 349.6 |
| LOX purge / propellant, worst point in burn [%] | 102.5 |
| LOX purge meets max propellant flow | yes |
| Fuel purge path CdA [mm2] | 1.75 |
| Fuel max propellant volumetric flow [L/s] | 0.555 |
| Fuel purge flow at that moment [L/s] | 0.457 |
| Fuel purge / propellant at peak [%] | 82.4 |
| Fuel feed pressure at peak [psia] | 384.9 |
| Fuel purge / propellant, worst point in burn [%] | 82.4 |
| Fuel purge meets max propellant flow | NO |
| Both purges at once vs regulator capacity [%] | 200.0 |
| Purge sequence | LOX Purge 2 s from T+0, Fuel Purge 2 s from T+0.1 |
| Purge requirement | >= 100 % of peak propellant volume flow |
| COPV when purges open [psia] | 1380.1 |
| COPV temperature when purges open [K] | 209.6 |
| LOX Purge required flow [L/s] at [psia] | 0.490 at 349.6 |
| Fuel Purge required flow [L/s] at [psia] | 0.555 at 384.9 |
| Minimum CV-4 / CV-5 CdA [mm2] | none up to 100 mm2: even at the highest set the COPV allows, a purge falls short. Usually the regulator's capacity (see regulator_sizing.csv / Hopper_Regulator_Sizing.py) |
| Regulator set needed at current CV-4 / CV-5 8.5 mm2 [psig] | 1272.0 |

# Free-Flight Sim: Inputs

## Thrust profile

| Input | Value | Notes |
| --- | --- | --- |
| Profile | 52m_free | GNC 6DOF reference, 52 m free flight (29.5 s) |
| Source | thrust_curves/thrust_profile_52m_free.csv |  |
| Interpolation | step | step: each point's flow held until the next |
| Throttle slew limit | 180 lbf/s | EDIT thrust slew limit (GNC fluid-network slides: 180 lbf/s sim, 90 lbf/s MC). null = step in step_smoothing_s (vehicle_sim-like) |
| Shortest throttle move | 0.01 s |  |
| CSV resample spacing | 0.1 s |  |

## Engine

| Input | Value | Notes |
| --- | --- | --- |
| Target mixture ratio | 1.25 | EDIT target O/F (user, 2026-10-30; was 1.8 from vehicle_sim) |
| Throat area | 1.042 in² |  |
| Exit area | 3.651 in² | Expansion ratio 3.50 |
| c* efficiency | 0.8 at 100 psia → 0.85 at 300 psia | Linear, flat outside |
| Cf efficiency | 0.94 at 100 psia → 0.97 at 250 psia | Linear, flat outside |
| Ambient pressure | 14.8 psia |  |
| Propellants / pressurant | LOX / RP-1 / Nitrogen |  |

## Start-up sequence

| Event | Time | Notes |
| --- | --- | --- |
| OIMV (torch ox) opens | T+0.08 s |  |
| FIMV (torch fuel) opens | T+0.19 s |  |
| FTV (main fuel) starts opening | T+0.3 s | Startup (user, 2026-11-10): FTV starts opening after the 100 ms torch Pc check; its flow ramps to the first profile point by main_ignition_s |
| OTV (main ox) starts opening | T+0.45 s | Startup (user, 2026-11-10): OTV starts opening (fuel lead) |
| Main ignition (profile start) | T+0.7 s | Startup (user, 2026-11-10): both flows reach the first thrust-profile point; profile time 0 |
| TFPV (torch purge) opens | T+0.73 s |  |
| OIMV / FIMV close | T+0.755 s |  |

## Shutdown sequence

Shutdown request = the end of the thrust profile.

| Event | Time | Notes |
| --- | --- | --- |
| OTV closes | shutdown + 0 s | Shutdown (user, 2026-11-10): OTV closes first at the shutdown request (end of the thrust profile), fuel-rich shutdown |
| FTV closes | shutdown + 0.1 s | Shutdown (user, 2026-11-10): FTV closes 0.1 s after OTV |
| OTV / FTV closing time | 0.2 s | GUESS OTV/FTV closing time |
| LOX Purge open | shutdown + 0 s for 2 s | OPV opens as the main oxidizer valve closes |
| Fuel Purge open | shutdown + 0.1 s for 2 s | FPV opens as the main fuel valve closes |
| TFPV closes | 0.1 s after the last main purge |  |
| Simulated after the last valve closes | 0.3 s |  |

## COPVs

| Bottle | Volume | Pressure | Temperature | Min temperature | Outlet CdA | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| COPV 1 | 9 L | 4503.2 psia | 293 K | 220 K | 100 mm² | GUESS outlet CdA and minimum temperature (Vedanth, 2026-11-05) |
| COPV 2 | 9 L | 4503.2 psia | 293 K | 220 K | 100 mm² | GUESS outlet CdA and minimum temperature (Vedanth, 2026-11-05) |

| Input | Value | Notes |
| --- | --- | --- |
| COPV model | lumped | lumped: all bottles as one volume (user, needs identical bottles). header: each bottle feeding a common header (Vedanth, 2026-11-05), for bottles that differ; stiff, so it needs maximum_step_s near 0.005 |
| Bang-bang sample period | 0.05 s | Bang-bang sampled once per time step, as a flight computer sampling its transducer |

## Tanks

| Tank | Propellant | Volume | Load | Temperature | Set pressure | Ullage gas | Collapse factor | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LOX Tank | LOX | 15.867 L | 15.998 kg | 90 K | 500 -5/+5 psia | 150 K | 2 | EDIT collapse factor assumed (user, 2026-10-27); Luca to supply conservative + Vespula values. Tank dimensions table (2026-10-30): 8.5 in OD, 10 % ullage; load = 45.54 - 10.27 lb wet - dry |
| Fuel Tank | RP-1 | 17.929 L | 13.000 kg | 293 K | 500 -5/+5 psia | 293 K | 1.1 | EDIT collapse factor assumed (user, 2026-10-27); Luca to supply conservative + Vespula values. Tank dimensions table (2026-10-30): 8.5 in OD, 10 % ullage; load = 42.22 - 13.56 lb wet - dry |

## LOX pressurant leg (COPV → LOX ullage)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| LOX Pressurant Line | Gas line | 0.25 in OD x 0.035 in wall (0.180 in ID), 48 in long | EDIT (Vedanth, 2026-11-05; was 4.92 mm ID, 1.0 m) |
| OBB | Bang-bang solenoid | 18.68 mm² | Instant solenoid (user assumption); set open/close times to model travel, which caps the time step at 1/10 of a stroke |
| OR-1 | Orifice | 7 mm² | EDIT (user, 2026-11-03; was 2.04 from vehicle_sim) |
| PROCK | Check valve | 12 mm² | EDIT |

## Fuel pressurant leg (COPV → fuel ullage)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| Fuel Pressurant Line | Gas line | 0.25 in OD x 0.035 in wall (0.180 in ID), 48 in long | EDIT (Vedanth, 2026-11-05; was 4.92 mm ID, 1.0 m) |
| FBB | Bang-bang solenoid | 18.68 mm² | Instant solenoid (user assumption); set open/close times to model travel, which caps the time step at 1/10 of a stroke |
| OR-2 | Orifice | 3 mm² | EDIT (user, 2026-11-03; was 1.02 from vehicle_sim) |
| PRFCK | Check valve | 12 mm² | EDIT |

## LOX feed

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| LOX Tank to OTV Line | Line | 0.375 in OD x 0.01 in wall (0.355 in ID), 48 in long | EDIT 3/8 in x 0.01 in wall (user, 2026-11-03). Length EDIT (Vedanth, 2026-11-05; was 1.0 m) |
| OTV | Throttle valve | 50 mm² full open | EDIT full-open CdA (user, 2026-11-03). Moves at the sequence thrust slew limit; set open/close times to add valve travel |
| OTV to LOX Injector Line | Line | 0.375 in OD x 0.01 in wall (0.355 in ID), 10 in long | EDIT 3/8 in x 0.01 in wall (user, 2026-11-03). Length EDIT (Vedanth, 2026-11-05; was 0.2032 m) |
| OINJ | Injector | 16.6 mm² | EDIT LOX stem (user, 2026-10-30; was 15.3 from vehicle_sim) |

## Fuel feed

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| Fuel Tank to FTV Line | Line | 0.375 in OD x 0.01 in wall (0.355 in ID), 12 in long | EDIT 3/8 in x 0.01 in wall (user, 2026-11-03; was 1/4 in from vehicle_sim) |
| FTV | Throttle valve | 50 mm² full open | EDIT full-open CdA (user, 2026-11-03). Moves at the sequence thrust slew limit; set open/close times to add valve travel |
| FTV to Regen Line | Line | 0.375 in OD x 0.01 in wall (0.355 in ID), 12 in long | EDIT 3/8 in x 0.01 in wall (user, 2026-11-03). Length EDIT (Vedanth, 2026-11-05; was 0.2032 m) |
| Regen Jacket | Regen | 30 mm²; +82 K at 150 psia, +82 K at 300 psia | EDIT CdA (user, 2026-11-03; was 25). EDIT heating: 293 K in, 375 K at the injector inlet (user, 2026-11-04); Vedanth's curve was +51 K at 150 psia, +44 K at 300 psia. Passage length is flow inertia, needed so the solver does not stall at shutdown |
| FINJ | Injector | 15 mm² | EDIT fuel manifold (user, 2026-10-30; was 10.6 from vehicle_sim) |

## Torch igniter

Enabled in the flight sim: yes. Valve strokes are each valve's open / close time.

### Torch oxidizer leg (OTV Inlet → Torch Chamber)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| OIMV | Valve | 0.68 mm² | Commands are in torch.standalone_purge.sequence. Propellant flow/combustion is excluded from the independent cold-gas purge run. |
| OICK | Check valve | 5 mm² | Oxidizer igniter check valve; INPUT CdA. |
| OIORF | Orifice | 8 mm² | User provisional CdA 8 mm2. |
| OTHX | Tube | 0.25 in OD x 0.065 in wall (0.120 in ID), 60 in long | Ox torch heat exchanger treated hydraulically as a tube. INPUT length and ID or OD/wall. Roughness and zero height are placeholders. No vaporization/phase-change model is implied. |
| Torch Ox Injector | Injector | 0.4398 mm² | INPUT torch oxidizer injector CdA. |

### Torch fuel leg (FTV Inlet → Torch Chamber)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| FIMV | Valve | 2.2 mm² | Commands are in torch.standalone_purge.sequence. The independent purge run plots travel but excludes propellant flow/combustion. |
| FICK | Check valve | 5 mm² | Fuel igniter check valve; INPUT CdA. |
| FIORF | Orifice | 8 mm² | User provisional CdA 8 mm2. Purge tees downstream, upstream of the fuel injector. |
| Torch Fuel Injector | Injector | 0.8051 mm² | INPUT torch fuel injector CdA. |

### Torch purge leg (PR-REG Outlet → FIORF Outlet)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| TFPV | Valve | 2.2 mm² | Renamed FIPV to TFPV. Commands are in torch.standalone_purge.sequence; open/close times here describe hardware travel. |
| FIPCK | Check valve | 5 mm² | Fuel purge igniter check valve; INPUT CdA. |
| FIPORF | Orifice | 50 mm² | Current user-entered CdA 50 mm2; earlier exploratory cases used 8 mm2. |

### Torch chamber and throat

| Input | Value | Notes |
| --- | --- | --- |
| Torch chamber volume | 0.0019 L | Solved as zero volume in the flight sim |
| Torch throat area | 0.006504 in² | Geometric throat area, no separate discharge coefficient. Independent purge run uses ideal bidirectional compressible throat flow. |

## Purge system

| Regulator | Source | Set | CdA | Supply effect | Droop | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| PR-REG | COPV | 750 psig | 0.8489 mm² | +17 psi per 1000 psi drop | 10 psi | Droop GUESS (not on the datasheet): outlet falls 10 psi from lock-up to fully open. Aqua Environment 1026 dome-loaded reducing regulator (user, 2026-11-10): Cv 0.05 (0.07 in orifice) = 0.849 mm2, 6000 psi rated, 0-6000 psi out, 1/4 in NPT, 17 psi outlet rise per 1000 psi inlet drop. Set 750 psig EDIT (applies at the starting COPV pressure); dome loading not modelled. Was 5 mm2 GUESS |

### ORV-1 Relief (PR-REG Outlet → ambient)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| ORV-1 | Relief valve | 2 mm², cracks 825 psig | EDIT GUESS: 10 % above PR-REG set |

### LOX Purge (PR-REG Outlet → OTV Outlet)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| LOX Purge Line | Gas line | 0.25 in OD x 0.035 in wall (0.180 in ID), 39.37 in long | EDIT 1/4 in x 0.035 in wall (user, 2026-11-08; was 0.028 in), 4.57 mm ID. Length GUESS |
| OPV | Valve | 1.86 mm² | EDIT candidate Chinese solenoid (user, 2026-10-31) |
| CV-5 | Check valve | 8.5 mm² | EDIT GUESS Cv ~0.5, typical 1/4 in poppet check (CdA = Cv/38 in^2) |

### Fuel Purge (PR-REG Outlet → FTV Outlet)

| Component | Type | Value | Notes |
| --- | --- | --- | --- |
| Fuel Purge Line | Gas line | 0.25 in OD x 0.035 in wall (0.180 in ID), 39.37 in long | EDIT 1/4 in x 0.035 in wall (user, 2026-11-08; was 0.028 in), 4.57 mm ID. Length GUESS |
| FPV | Valve | 1.86 mm² | EDIT candidate Chinese solenoid (user, 2026-10-31) |
| CV-4 | Check valve | 8.5 mm² | EDIT GUESS Cv ~0.5, typical 1/4 in poppet check (CdA = Cv/38 in^2) |

## Vents, reliefs and drain

| Path | From → to | Component | Type | Value | Notes |
| --- | --- | --- | --- | --- | --- |
| PRVNT Vent | COPV → ambient | PRVNT | Valve | 1.07 mm² | Marotta PLV34 (info); 1x PLV34 = 1.07 mm^2 (sheet FRV tab) |
| PRRV Relief | COPV → ambient | PRRV | Relief valve | 5 mm², cracks 4950 psig | EDIT GUESS: Swagelok SS-4R3A (info), sheet tab empty |
| OVNT Vent | LOX Ullage → ambient | OVNT | Valve | 9.4 mm² | Jaksa 2/2 cryo solenoid (info), 9.4 mm^2 (sheet OVNT) |
| ORV Relief | LOX Ullage → ambient | ORV | Relief valve | 30.8 mm², cracks 550 psig | RegO PRV9432T550 (info); target CdA (sheet ORV) |
| OFD Drain | LOX Tank Liquid → ambient | OFD | Valve | 9.4 mm² | Same valve as OVNT (Jaksa 2/2 cryo solenoid) |
| FVNT Vent | Fuel Ullage → ambient | FVNT | Valve | 15 mm² | UNOX 3/8 NO solenoid (info); target CdA (sheet FVNT) |
| FRV Relief | Fuel Ullage → ambient | FRV | Relief valve | 14.67 mm², cracks 550 psig | EDIT PRV9432T550 = 550 psig (info); sheet FRV says 525 |

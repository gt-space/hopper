| Path | Tag | Type | CdA [mm2] | Details | Source / note |
|---|---|---|---|---|---|
| Pressurant | COPV | COPV |  | 3 x 9 L, 4503 psia, 293 K Nitrogen |  |
| LOX Tank | LOX Tank | Tank |  | 20 L, 19.43 kg LOX, set 500 psia +5/-5, collapse factor 2 |  |
| Fuel Tank | Fuel Tank | Tank |  | 16 L, 9.7 kg RP-1, set 500 psia +5/-5, collapse factor 1.1 |  |
| LOX pressurant leg | LOX Pressurant Line | Gas line | 7.55 | ID 4.92 mm, L 1.000 m; lumped (fully rough) | meter, EDIT |
| LOX pressurant leg | OBB | Bang-bang solenoid | 18.68 | opened/closed by the tank bang-bang band |  |
| LOX pressurant leg | OR-1 | Orifice | 2.04 |  | EDIT |
| LOX pressurant leg | PROCK | Check valve | 12.00 |  | EDIT |
| Fuel pressurant leg | Fuel Pressurant Line | Gas line | 7.55 | ID 4.92 mm, L 1.000 m; lumped (fully rough) | meter, EDIT |
| Fuel pressurant leg | FBB | Bang-bang solenoid | 18.68 | opened/closed by the tank bang-bang band |  |
| Fuel pressurant leg | OR-2 | Orifice | 1.02 |  | EDIT |
| Fuel pressurant leg | PRFCK | Check valve | 12.00 |  | EDIT |
| LOX feed | LOX Tank to OTV Line | Line |  | ID 8.10 mm, L 1.000 m |  |
| LOX feed | OTV | Throttle valve | 46.41 | max CdA; CdA set to hit the flow target | EDIT |
| LOX feed | OTV to LOX Injector Line | Line |  | ID 8.10 mm, L 0.203 m |  |
| LOX feed | OINJ | Injector | 15.30 | passage length 0.010 m |  |
| Fuel feed | Fuel Tank to FTV Line | Line |  | ID 4.93 mm, L 0.305 m |  |
| Fuel feed | FTV | Throttle valve | 17.16 | max CdA; CdA set to hit the flow target | EDIT: 90 % of tube area (17.2 mm^2) |
| Fuel feed | FTV to Regen Line | Line |  | ID 4.93 mm, L 0.203 m |  |
| Fuel feed | Regen Jacket | Orifice | 25.00 | passage length 0.300 m | EDIT |
| Fuel feed | FINJ | Injector | 10.60 | passage length 0.010 m |  |
| COPV side path | PR-REG | Regulator | 5.00 | ideal, set 500 psig | EDIT GUESS: not in sheet |
| PRVNT Vent: COPV -> ambient | PRVNT | Valve | 1.07 | closed all run | Marotta PLV34 (info); 1x PLV34 = 1.07 mm^2 (sheet FRV tab) |
| PRRV Relief: COPV -> ambient | PRRV | Relief valve | 5.00 | crack 4950 psig, full open 5445 psig | EDIT GUESS: Swagelok SS-4R3A (info), sheet tab empty |
| ORV-1 Relief: PR-REG Outlet -> ambient | ORV-1 | Relief valve | 2.00 | crack 600 psig, full open 660 psig | EDIT GUESS |
| LOX Purge: PR-REG Outlet -> OTV Outlet | OPV | Valve | 2.00 | closed all run | EDIT GUESS CdA |
| LOX Purge: PR-REG Outlet -> OTV Outlet | CV-5 | Check valve | 5.00 |  | EDIT GUESS CdA |
| Fuel Purge: PR-REG Outlet -> FTV Outlet | FPV | Valve | 2.00 | closed all run | EDIT GUESS CdA |
| Fuel Purge: PR-REG Outlet -> FTV Outlet | CV-4 | Check valve | 5.00 |  | EDIT GUESS CdA |
| Torch Fuel Purge: PR-REG Outlet -> Main Chamber | TFPV | Valve | 2.00 | closed all run | EDIT GUESS CdA |
| Torch Fuel Purge: PR-REG Outlet -> Main Chamber | CV-3 | Check valve | 5.00 |  | EDIT GUESS CdA |
| OVNT Vent: LOX Ullage -> ambient | OVNT | Valve | 9.40 | closed all run | Jaksa 2/2 cryo solenoid (info), 9.4 mm^2 (sheet OVNT) |
| ORV Relief: LOX Ullage -> ambient | ORV | Relief valve | 30.80 | crack 550 psig, full open 605 psig | RegO PRV9432T550 (info); target CdA (sheet ORV) |
| OFD Drain: LOX Tank Liquid -> ambient | OFD | Valve | 9.40 | closed all run | same valve as OVNT (Jaksa 2/2 cryo solenoid) |
| FVNT Vent: Fuel Ullage -> ambient | FVNT | Valve | 15.00 | closed all run | UNOX 3/8 NO solenoid (info); target CdA (sheet FVNT) |
| FRV Relief: Fuel Ullage -> ambient | FRV | Relief valve | 14.67 | crack 550 psig, full open 605 psig | EDIT PRV9432T550 = 550 psig (info); sheet FRV says 525 |
| Not in the flow model | PRQD |  |  |  | Pressurant fill quick disconnect (disconnected in flight) |
| Not in the flow model | OQD |  |  |  | LOX fill quick disconnect below OFD (disconnected in flight) |
| Not in the flow model | FQD |  |  |  | Fuel fill quick disconnect, McMaster 52495K63 (disconnected in flight) |
| Not in the flow model | ODiff |  |  |  | Oxygen diffuser on the OVNT exhaust (sizing) |
| Not in the flow model | Torch spark igniter |  |  |  | Spark / exciter; torch combustion itself is not modelled (its propellant goes into the main chamber) |
| Not in the flow model | PRPT-1 and sensors |  |  |  | Instrumentation |
| Not in the flow model | OIMV / FIMV / HX-1 |  |  |  | Torch igniter feeds, switched off (TORCH_ENABLED = False) |

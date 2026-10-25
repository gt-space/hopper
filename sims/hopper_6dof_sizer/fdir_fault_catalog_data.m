function [severity, tier] = fdir_fault_catalog_data()
% FDIR_FAULT_CATALOG_DATA  Severity and response tier for all 102 faults.
%
% GENERATED FILE - do not edit by hand.
% Run make_fdir_catalog to rebuild it from fdir_fault_catalog.csv,
% which is the "Fault List" sheet of FaultConditions.xlsx.
%
% Severity 0 means the sheet lists the fault as N/A (not on the vehicle).
% Tier is severity - 1: 0 -> R0, 1 -> R1, 2 -> R2, 3 -> R3, 4 -> R4.
%
%#codegen

severity = [ ...
    3, ...   % F-01 GNC / Sensors: IMU loss (no data / stale data) -> R2
    3, ...   % F-02 GNC / Sensors: IMU drift or bias beyond estimate error limits (+/-1 deg attitude) -> R2
    2, ...   % F-03 GNC / Sensors: Magnetometer loss or interference -> R1
    3, ...   % F-04 GNC / Sensors: Barometer loss or plume/pressure-induced error -> R2
    3, ...   % F-05 GNC / Sensors: GPS/GNSS dropout -> R2
    2, ...   % F-06 GNC / Sensors: GNSS multipath / degraded position solution -> R1
    3, ...   % F-07 GNC / Sensors: LiDAR altimeter dropout (dust, reflective surfaces, plume) -> R2
    3, ...   % F-08 GNC / Sensors: Loss of one or more of 3+ LiDAR modules / LiDAR disagreement -> R2
    4, ...   % F-09 GNC / Sensors: State estimate (EKF) divergence -> R3
    3, ...   % F-10 GNC / Sensors: Position estimate error > +/-1 m / velocity error > +/-1 m/s -> R2
    3, ...   % F-11 GNC / Sensors: Vehicle mass estimate error > +/-1 kg -> R2
    1, ...   % F-12 GNC / Sensors: Pad load cell inaccuracy / failure prior to takeoff -> R0
    5, ...   % F-13 GNC / Control: Controller instability (poor tuning, positive feedback loop) -> R4
    4, ...   % F-14 GNC / Control: Limit cycling / "buzzing" from actuator deadband -> R3
    5, ...   % F-15 GNC / Control: Loss of control authority (TWR < 1) -> R4
    4, ...   % F-16 GNC / Control: Position error > 5 m -> R3
    4, ...   % F-17 GNC / Control: Attitude error > 5 deg / tilt cone exceedance (> 7 deg nominal, > 2 deg landing) -> R3
    4, ...   % F-18 GNC / Control: Body rate error > 5 deg/s -> R3
    3, ...   % F-19 GNC / Control: Velocity exceedance (> 8 m/s) -> R2
    4, ...   % F-20 GNC / Control: Excessive landing normal (> 1.5 m/s) or lateral (> 0.5 m/s) velocity -> R3
    4, ...   % F-21 GNC / Control: Landing outside 1 m radius / on pad hardware (LFS) -> R3
    4, ...   % F-22 GNC / Control: Wind disturbance exceeding 10 m/s design limit -> R3
    5, ...   % F-23 GNC / Guidance: Infeasible trajectory generated (G-FOLD failure) -> R4
    3, ...   % F-24 GNC / Guidance: Trajectory causes actuator saturation -> R2
    2, ...   % F-25 GNC / Dynamics: Slosh-induced oscillations (ox or fuel) -> R1
    2, ...   % F-26 GNC / Dynamics: Structural bending / resonance coupling with TVC -> R1
    2, ...   % F-27 GNC / Dynamics: Tail-wag instability from gimbal-induced CG shift -> R1
    2, ...   % F-28 GNC / Dynamics: Thrust misalignment (mounting error) -> R1
    5, ...   % F-29 TVC: TVC actuator failure (stuck / frozen) -> R4
    2, ...   % F-30 TVC: TVC actuator hard-over to +/-8 deg limit -> R1
    2, ...   % F-31 TVC: TVC rate saturation (< 15 deg/s achieved) -> R1
    4, ...   % F-32 TVC: TVC position error / encoder failure -> R3
    2, ...   % F-33 TVC: TVC latency exceeding 100 Hz loop requirement -> R1
    4, ...   % F-34 RCS: RCS fan (EDF) failure / loss of roll control -> R3
    3, ...   % F-35 RCS: RCS fan thrust saturation (> 16 N demand) -> R2
    5, ...   % F-36 Avionics: Flight computer crash / reset / watchdog timeout -> R4
    5, ...   % F-37 Avionics: Compute latency causing phase lag (loop < 100 Hz) -> R4
    5, ...   % F-38 Avionics: Engine controller (EC) failure -> R4
    5, ...   % F-39 Avionics: Ethernet link loss between vehicle boards -> R4
    3, ...   % F-40 Avionics: Loss of telemetry / comms with ground (FTEL/GTEL) -> R2
    5, ...   % F-41 Avionics: Loss of main avionics power (8S1P battery) -> R4
    5, ...   % F-42 Avionics: Loss of motor battery pack (throttle valve / TVC motors) -> R4
    4, ...   % F-43 Avionics: BMS fault (over/under voltage, overcurrent, overtemp) -> R3
    3, ...   % F-44 Avionics: Avionics enclosure overtemperature -> R2
    4, ...   % F-45 Avionics: Vibration-induced connector/harness failure -> R3
    5, ...   % F-46 Avionics: Umbilical fails to disconnect at launch -> R4
    2, ...   % F-47 Avionics: Premature umbilical disconnect before launch -> R1
    2, ...   % F-48 Avionics: Blackbox / SD card data logging failure -> R1
    2, ...   % F-49 Fluids / Sensors: Pressure transducer failure (OPT, FPT, PRPT, throttle valve PTs) -> R1
    1, ...   % F-50 Fluids / Sensors: Temperature sensor failure (OUTC, OTC, FTC, PRRTD) -> R0
    1, ...   % F-51 Fluids / Sensors: Sensor reading disagreement between redundant PTs -> R0
    2, ...   % F-52 Fluids: Ox tank over-pressure (> 500 psi MEOP) -> R1
    2, ...   % F-53 Fluids: Fuel tank over-pressure (> 500 psi MEOP) -> R1
    2, ...   % F-54 Fluids: Ox or fuel tank under-pressure / pressure collapse -> R1
    2, ...   % F-55 Fluids: Tank pressure oscillation > +/-30 psi -> R1
    3, ...   % F-56 Fluids: Bang-bang valve (OBB-1/OBB-2/FBB) stuck closed -> R2
    3, ...   % F-57 Fluids: Bang-bang valve stuck open -> R2
    3, ...   % F-58 Fluids: COPV over-pressure (> 4500 psi) -> R2
    5, ...   % F-59 Fluids: COPV pressure depletion before end of flight -> R4
    2, ...   % F-60 Fluids: COPV below minimum temperature (< 233 K) -> R1
    3, ...   % F-61 Fluids: Pressurant vent (PRVNT) failure to open/close -> R2
    3, ...   % F-62 Fluids: Ox vent (OVNT) failure to open/close -> R2
    3, ...   % F-63 Fluids: Fuel vent (FVNT) failure to open/close -> R2
    3, ...   % F-64 Fluids: Relief valve (ORV/FRV/PRRV) fails to relieve or leaks -> R2
    2, ...   % F-65 Fluids: Check valve (OCK/FCK) backflow / failure -> R1
    2, ...   % F-66 Fluids: Propellant or pressurant leak -> R1
    2, ...   % F-67 Fluids: Chinese solenoid (FIMV/OIMV/FVNT) overheating or failure -> R1
    2, ...   % F-68 Fluids: LOX boil-off / geysering in feed lines -> R1
    5, ...   % F-69 Fluids: Propellant depletion before landing -> R4
    2, ...   % F-70 Fluids: QD (OQD/FQD/PRQD/OIQD) leak or failure to release -> R1
    5, ...   % F-71 Engine: Oxidizer throttle valve (OTV) stuck / motor failure -> R4
    5, ...   % F-72 Engine: Fuel throttle valve (FTV) stuck / motor failure -> R4
    0, ...   % F-73 Engine: Throttle valve limit switch failure -> R0
    5, ...   % F-74 Engine: Throttle valve cryogenic binding / seat failure -> R4
    2, ...   % F-75 Engine: Throttle response too slow / pressure loop tracking error -> R1
    4, ...   % F-76 Engine: Mixture ratio excursion from 1.25 -> R3
    2, ...   % F-77 Engine: Chamber pressure off-nominal (vs. 300 psia) -> R1
    2, ...   % F-78 Engine: Thrust below required TWR (low throttle / underperformance) -> R1
    3, ...   % F-79 Engine: Throttle outside 40-110% range -> R2
    2, ...   % F-80 Engine: Low injector stiffness (< 15%) / combustion instability / chugging -> R1
    3, ...   % F-81 Engine: Hard start -> R2
    1, ...   % F-82 Engine: Ignition failure (torch igniter no-light) -> R0
    5, ...   % F-83 Engine: Flameout during flight -> R4
    3, ...   % F-84 Engine: Regen chamber burn-through / overheating -> R2
    1, ...   % F-85 Engine: Thrust cut failure at 0.05 m altitude -> R0
    1, ...   % F-86 Engine: Failed shutdown (propellant flow not stopped) -> R0
    2, ...   % F-87 Engine: Igniter GOX system failure (OIRV) -> R1
    2, ...   % F-88 Engine: OIVNT Failure -> R1
    3, ...   % F-89 Engine: OIREG Failure -> R2
    4, ...   % F-90 Structures: Landing leg failure / collapse -> R3
    5, ...   % F-91 Structures: Vehicle tip-over on landing -> R4
    5, ...   % F-92 Structures: Thrust structure failure -> R4
    3, ...   % F-93 Structures: Baffle failure / detachment -> R2
    4, ...   % F-94 GSE: Tether failure / excessive tether loading -> R3
    4, ...   % F-95 GSE: Tether snag / entanglement during hop -> R3
    3, ...   % F-96 GSE: Ground computer / GSAM failure -> R2
    3, ...   % F-97 GSE: Flame deflector / pad damage -> R2
    5, ...   % F-98 System: Loss of power (abort scenario) -> R4
    5, ...   % F-99 System: Loss of comms (abort scenario) -> R4
    5, ...   % F-100 System: Loss of control (abort scenario) -> R4
    5, ...   % F-101 System: Kill switch / abort command failure -> R4
    5, ...   % F-102 System: Flight outside approved range / airspace boundary -> R4
    ];

tier = [ ...
    2, 2, 1, 2, 2, 1, 2, 2, 3, 2, ...   % F-01..F-10
    2, 0, 4, 3, 4, 3, 3, 3, 2, 3, ...   % F-11..F-20
    3, 3, 4, 2, 1, 1, 1, 1, 4, 1, ...   % F-21..F-30
    1, 3, 1, 3, 2, 4, 4, 4, 4, 2, ...   % F-31..F-40
    4, 4, 3, 2, 3, 4, 1, 1, 1, 0, ...   % F-41..F-50
    0, 1, 1, 1, 1, 2, 2, 2, 4, 1, ...   % F-51..F-60
    2, 2, 2, 2, 1, 1, 1, 1, 4, 1, ...   % F-61..F-70
    4, 4, 0, 4, 1, 3, 1, 1, 2, 1, ...   % F-71..F-80
    2, 0, 4, 2, 0, 0, 1, 1, 2, 3, ...   % F-81..F-90
    4, 4, 2, 3, 3, 2, 2, 4, 4, 4, ...   % F-91..F-100
    4, 4 ...   % F-101..F-102
    ];
end


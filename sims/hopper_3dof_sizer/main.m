% File Additions
addpath('./sizing')
addpath('./inputs')
addpath('./propulsion')
addpath('./dynamics')

% Main Vehicle Sizing Section
% --- Tank, AVI, and Engine Sized prior to primary structure calculations
IN = mission_inputs();

TANKS = size_tanks(IN);
AVI = size_avionics(IN);
ENGINE = size_engine(IN);

VEH.mass.tanks = TANKS.singular.total_mass;
VEH.mass.ox_tank = TANKS.singular.oxidizer.mass;
VEH.mass.fu_tank = TANKS.singular.fuel.mass;
VEH.h.ox_tank = TANKS.singular.oxidizer.h;
VEH.h.fu_tank = TANKS.singular.fuel.h;
VEH.mass.engine = ENGINE.mass;
VEH.mass.avi = AVI.mass;
VEH.mass.dry = VEH.mass.engine + VEH.mass.tanks + VEH.mass.avi + IN.structures.payload_mass + IN.structures.extra_payload_mass;
% PROP = size_propellant(IN, VEH, ENGINE); need engine code to rum this

VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass ;

STRUCT = size_structures(IN, VEH);

% VEHICLE SIZING OUTPUTS
OUT = Outputs(IN, VEH, TANKS, STRUCT);
fprintf('\n===== VEHICLE SUMMARY =====\n');
fprintf('Dry Mass: %.2f kg\n', OUT.Vehicle.DryMass);
fprintf('Wet Mass: %.2f kg\n', OUT.Vehicle.WetMass);
fprintf('Propellant Mass Percentage: %.2f\n', OUT.Vehicle.MassRatio)
fprintf('Initial TWR: %.2f\n', OUT.Vehicle.InitialTWR);
fprintf('Final TWR: %.2f\n', OUT.Vehicle.FinalTWR);

fprintf('\n=== Propellant ===\n');
fprintf('Oxidizer Mass: %.2f kg\n', IN.propulsion.oxidizer_mass);
fprintf('Fuel Mass: %.2f kg\n', IN.propulsion.fuel_mass);

fprintf('\n=== Tanks ===\n');
fprintf('Ox Tank Volume: %.4f L\n', OUT.Structures.OxTankVolume * 1000);
fprintf('Fuel Tank Volume: %.4f L\n', OUT.Structures.FuelTankVolume * 1000);
%fprintf('OX Tank Mass: %.2f kg\n', OUT.Structures.FuelTankVolume * 1000);
%fprintf('FU Tank Mass: %.2f kg\n', OUT.Structures.FuelTankVolume * 1000);

fprintf('\n=== Structures ===\n');
fprintf('Landing Legs + Intertank: %.2f kg\n', ...
    OUT.Vehicle.MassDistribution.Structures);

fprintf('\n=== Avionics ===\n');
fprintf('Battery Capacity: %.4f Ah\n', AVI.battery_capacity);
fprintf('# of Cells: %.4f \n', AVI.num_cells);
fprintf('Mass: %.4f kg\n', AVI.mass);

fprintf('\n=== Engines ===\n');
fprintf('Injector Mass: %.2f kg\n', ENGINE.inj_mass);
fprintf('TCA Mass: %.2f kg\n', ENGINE.TCA_mass);
fprintf('TVC Mass: %.2f kg\n', ENGINE.TVC_mass);
fprintf('Engine Mass: %.2f kg\n', ENGINE.mass);

fprintf('\n=== Payload ===\n');
fprintf('CPLC Payload Mass: %.2f kg\n', IN.structures.payload_mass);
fprintf('Remaining Payload Mass: %.2f kg\n', IN.structures.extra_payload_mass);

% File Additions
addpath('./sizing')
addpath('./inputs')
addpath('./propulsion')
addpath('./dynamics')

% Main Vehicle Sizing Section
IN = mission_inputs();

TANKS = size_tanks(IN);
AVI = size_avionics(IN, vehicle);

VEH.mass.tanks = TANKS.singular.total_mass;
VEH.mass.engine = 12; % placeholder
VEH.mass.avi = AVI.mass;
VEH.mass.dry = VEH.mass.engine + VEH.mass.tanks + IN.structures.payload_mass + VEH.mass.avi;
% PROP = size_propellant(IN, VEH, ENGINE); need engine code to rum this

VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass ;

STRUCT = size_structures(IN, VEH);

OUT = Outputs(IN, VEH, TANKS, STRUCT);
fprintf('\n===== VEHICLE SUMMARY =====\n');

fprintf('Dry Mass: %.2f kg\n', OUT.Vehicle.DryMass);
fprintf('Wet Mass: %.2f kg\n', OUT.Vehicle.WetMass);

fprintf('Initial TWR: %.2f\n', OUT.Vehicle.InitialTWR);
fprintf('Final TWR: %.2f\n', OUT.Vehicle.FinalTWR);

fprintf('\n--- Propellant ---\n');
fprintf('Oxidizer Mass: %.2f kg\n', IN.propulsion.oxidizer_mass);
fprintf('Fuel Mass: %.2f kg\n', IN.propulsion.fuel_mass);

fprintf('\n--- Tank Volumes ---\n');
fprintf('Ox Tank Volume: %.4f m^3\n', OUT.Structures.OxTankVolume);
fprintf('Fuel Tank Volume: %.4f m^3\n', OUT.Structures.FuelTankVolume);

fprintf('\n--- Structures ---\n');
fprintf('Landing Legs + Intertank: %.2f kg\n', ...
    OUT.Vehicle.MassDistribution.Structures);

fprintf('\n--- Avionics ---\n');
fprintf('Battery Capacity: %.4f Ah\n', AVI.battery_capacity);
fprintf('# of Cells: %.4f \n', AVI.num_cells);
fprintf('Mass: %.4f kg\n', AVI.mass);


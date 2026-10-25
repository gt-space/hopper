
% Main Vehicle Sizing Section
% --- Tank, AVI, and Engine Sized prior to primary structure calculations
IN = mission_inputs();

TANKS = size_tanks(IN);
AVI = size_avionics(IN);
ENGINE = size_engine(IN);

% Dry mass override: the vehicle mass budget until the CAD mass-properties
% export exists. Includes the legs and intertank.
dry_mass_override = 116.86;
VEH.mass.dry = dry_mass_override;

VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass;

[STRUCT, VEH] = size_structures(IN, VEH);
VEH.mass.dry = dry_mass_override;   % size_structures adds legs/intertank; already in the override

r_leg = [
    0 IN.legs.feet_radial_dist 0
    0 -IN.legs.feet_radial_dist 0
    0 0 IN.legs.feet_radial_dist
    0 0 -IN.legs.feet_radial_dist
]';
r_leg_1 = r_leg(:,1)';   % per-leg names still read by the v2/v3 models
r_leg_2 = r_leg(:,2)';
r_leg_3 = r_leg(:,3)';
r_leg_4 = r_leg(:,4)';

VEH.k_engine = (IN.engine_cg/IN.mission.max_vehicle_height)^2;


% VEHICLE SIZING OUTPUTS
OUT = Outputs(IN, VEH, TANKS, STRUCT);
fprintf('\n===== VEHICLE SUMMARY =====\n');
fprintf('Dry Mass: %.2f kg\n', OUT.Vehicle.DryMass);
fprintf('Wet Mass: %.2f kg\n', OUT.Vehicle.WetMass);
fprintf('Propellant Mass Percentage: %.2f\n', OUT.Vehicle.MassRatio)
fprintf('Initial TWR: %.2f\n', OUT.Vehicle.InitialTWR);
fprintf('Final TWR: %.2f\n', OUT.Vehicle.FinalTWR);
fprintf('Min Thrust: %.3f N\n', OUT.Vehicle.minThrust)
fprintf('Max Thrust: %.3f N\n', OUT.Vehicle.maxThrust)

fprintf('\n=== Propellant ===\n');
fprintf('Oxidizer Mass: %.2f kg\n', IN.propulsion.oxidizer_mass);
fprintf('Fuel Mass: %.2f kg\n', IN.propulsion.fuel_mass);
fprintf('Oxidizer Tank Pressure: %.2f psia\n', IN.propulsion.ox_press / 6894.76);
fprintf('Fuel Tank Pressure: %.2f psia\n', IN.propulsion.fu_press / 6894.76);

fprintf('\n=== Tanks ===\n');
fprintf('Ox Tank Volume: %.4f L\n', OUT.Structures.OxTankVolume * 1000);
fprintf('Fuel Tank Volume: %.4f L\n', OUT.Structures.FuelTankVolume * 1000);
fprintf('Ox Tank Mass: %.4f kg\n',TANKS.singular.oxidizer.mass);
fprintf('Ox Tank Height: %.4f m\n',TANKS.singular.oxidizer.h);
fprintf('Fuel Tank Mass: %.4f kg\n',TANKS.singular.fuel.mass);
fprintf('Fuel Tank Height: %.4f m\n',TANKS.singular.fuel.h);
fprintf('Total Tank Mass: %.4f kg\n',TANKS.singular.total_mass);
fprintf('Tank Radius: %.4f m\n',TANKS.singular.radius);

fprintf('\n=== COPV ===\n');
fprintf('COPV Mass: %.2f kg\n', IN.press.copv.mass);
fprintf('COPV Total Gas Mass: %.2f kg\n', IN.press.copv.gas_mass);
fprintf('COPV Amount: %d\n', IN.press.copv.amount)
fprintf('COPV Max Volume: %.3f m^3\n', IN.press.copv.max_volume);
fprintf('COPV Max Pressure: %.3f psia \n', IN.press.copv.max_pressure / 6894.76);

fprintf('\n=== Structures ===\n');
fprintf('Landing Legs + Intertank: %.2f kg\n', ...
    OUT.Vehicle.MassDistribution.Structures);

fprintf('\n=== Avionics ===\n');
fprintf('Battery Capacity: %.4f Ah\n', AVI.battery_capacity);
fprintf('# of Cells: %.4f \n', AVI.num_cells);
fprintf('Mass: %.4f kg\n', AVI.mass);

fprintf('\n=== Engines ===\n');
fprintf('C Star Efficiency: %.2f \n', IN.propulsion.eta_cstar);
fprintf('Throat Area: %.3f in^2\n', IN.propulsion.At / 0.00064516);
fprintf('Expansion Ratio: %.3f \n', IN.propulsion.eps);
fprintf('Injector Mass: %.2f kg\n', ENGINE.inj_mass);
fprintf('TCA Mass: %.2f kg\n', ENGINE.TCA_mass);
fprintf('TVC Mass: %.2f kg\n', ENGINE.TVC_mass);
fprintf('Engine Mass: %.2f kg\n', ENGINE.mass);

fprintf('\n=== Payload ===\n');
fprintf('CPLC Payload Mass: %.2f kg\n', IN.structures.payload_mass);
% fprintf('Remaining Payload Mass: %.2f kg\n', IN.structures.extra_payload_mass);

function result = mc_main(scenario, sim_out)

% File Additions
addpath('./sizing')
addpath('./inputs')
addpath('./propulsion')
addpath('./dynamics')

% Main Vehicle Sizing Section
% --- Tank, AVI, and Engine Sized prior to primary structure calculations
IN = mission_inputs(scenario);

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
VEH.mass.press = IN.press.copv.mass * IN.press.copv.amount + IN.press.copv.gas_mass;
VEH.mass.dry = IN.mass_factor * VEH.mass.engine + VEH.mass.tanks + VEH.mass.avi + IN.structures.payload_mass + IN.structures.extra_payload_mass + VEH.mass.press;
% PROP = size_propellant(IN, VEH, ENGINE); need engine code to rum this

VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass ;

[STRUCT, VEH] = size_structures(IN, VEH);

r_leg_1 = [0 IN.legs.feet_radial_dist 0];
r_leg_2 = [0 -IN.legs.feet_radial_dist 0];
r_leg_3 = [0 0 IN.legs.feet_radial_dist];
r_leg_4 = [0 0 -IN.legs.feet_radial_dist];

VEH.k_engine = (IN.engine_cg/1.62)^2;

% VEHICLE SIZING OUTPUTS
OUT = Outputs(IN, VEH, TANKS, STRUCT);
% fprintf('\n===== VEHICLE SUMMARY =====\n');
% fprintf('Dry Mass: %.2f kg\n', OUT.Vehicle.DryMass);
% fprintf('Wet Mass: %.2f kg\n', OUT.Vehicle.WetMass);
% fprintf('Propellant Mass Percentage: %.2f\n', OUT.Vehicle.MassRatio)
% fprintf('Initial TWR: %.2f\n', OUT.Vehicle.InitialTWR);
% fprintf('Final TWR: %.2f\n', OUT.Vehicle.FinalTWR);
% fprintf('Min Thrust: %.3f N\n', OUT.Vehicle.minThrust)
% fprintf('Max Thrust: %.3f N\n', OUT.Vehicle.maxThrust)
% 
% fprintf('\n=== Propellant ===\n');
% fprintf('Oxidizer Mass: %.2f kg\n', IN.propulsion.oxidizer_mass);
% fprintf('Fuel Mass: %.2f kg\n', IN.propulsion.fuel_mass);
% fprintf('Oxidizer Tank Pressure: %.2f psia\n', IN.propulsion.ox_press / 6894.76);
% fprintf('Fuel Tank Pressure: %.2f psia\n', IN.propulsion.fu_press / 6894.76);
% 
% fprintf('\n=== Tanks ===\n');
% fprintf('Ox Tank Volume: %.4f L\n', OUT.Structures.OxTankVolume * 1000);
% fprintf('Fuel Tank Volume: %.4f L\n', OUT.Structures.FuelTankVolume * 1000);
% fprintf('Ox Tank Mass: %.4f kg\n',TANKS.singular.oxidizer.mass);
% fprintf('Ox Tank Height: %.4f m\n',TANKS.singular.oxidizer.h);
% fprintf('Fuel Tank Mass: %.4f kg\n',TANKS.singular.fuel.mass);
% fprintf('Fuel Tank Height: %.4f m\n',TANKS.singular.fuel.h);
% fprintf('Total Tank Mass: %.4f kg\n',TANKS.singular.total_mass);
% fprintf('Tank Radius: %.4f m\n',TANKS.singular.radius);
% 
% fprintf('\n=== COPV ===\n');
% fprintf('COPV Mass: %.2f kg\n', IN.press.copv.mass);
% fprintf('COPV Total Gas Mass: %.2f kg\n', IN.press.copv.gas_mass);
% fprintf('COPV Amount: %d\n', IN.press.copv.amount)
% fprintf('COPV Max Volume: %.3f m^3\n', IN.press.copv.max_volume);
% fprintf('COPV Max Pressure: %.3f psia \n', IN.press.copv.max_pressure / 6894.76);
% 
% fprintf('\n=== Structures ===\n');
% fprintf('Landing Legs + Intertank: %.2f kg\n', ...
%     OUT.Vehicle.MassDistribution.Structures);
% 
% fprintf('\n=== Avionics ===\n');
% fprintf('Battery Capacity: %.4f Ah\n', AVI.battery_capacity);
% fprintf('# of Cells: %.4f \n', AVI.num_cells);
% fprintf('Mass: %.4f kg\n', AVI.mass);
% 
% fprintf('\n=== Engines ===\n');
% fprintf('C Star Efficiency: %.2f \n', IN.propulsion.eta_cstar);
% fprintf('Throat Area: %.3f in^2\n', IN.propulsion.At / 0.00064516);
% fprintf('Expansion Ratio: %.3f \n', IN.propulsion.eps);
% fprintf('Injector Mass: %.2f kg\n', ENGINE.inj_mass);
% fprintf('TCA Mass: %.2f kg\n', ENGINE.TCA_mass);
% fprintf('TVC Mass: %.2f kg\n', ENGINE.TVC_mass);
% fprintf('Engine Mass: %.2f kg\n', ENGINE.mass);
% 
% fprintf('\n=== Payload ===\n');
% fprintf('CPLC Payload Mass: %.2f kg\n', IN.structures.payload_mass);
% fprintf('Remaining Payload Mass: %.2f kg\n', IN.structures.extra_payload_mass);


%OUT = Outputs(IN, VEH, TANKS, STRUCT);

result.scenario = scenario;

result.vehicle.dry_mass    = OUT.Vehicle.DryMass;
result.vehicle.wet_mass    = OUT.Vehicle.WetMass;
result.vehicle.mass_ratio  = OUT.Vehicle.MassRatio;
result.vehicle.twr_initial = OUT.Vehicle.InitialTWR;
result.vehicle.twr_final   = OUT.Vehicle.FinalTWR;
result.vehicle.thrust_min  = OUT.Vehicle.minThrust;
result.vehicle.thrust_max  = OUT.Vehicle.maxThrust;

result.propellant.ox_mass      = IN.propulsion.oxidizer_mass;
result.propellant.fuel_mass    = IN.propulsion.fuel_mass;
result.propellant.ox_press_psi = IN.propulsion.ox_press / 6894.76;
result.propellant.fu_press_psi = IN.propulsion.fu_press / 6894.76;

result.tanks.ox_volume_L = OUT.Structures.OxTankVolume * 1000;
result.tanks.fu_volume_L = OUT.Structures.FuelTankVolume * 1000;
result.tanks.ox_mass     = TANKS.singular.oxidizer.mass;
result.tanks.ox_height   = TANKS.singular.oxidizer.h;
result.tanks.fu_mass     = TANKS.singular.fuel.mass;
result.tanks.fu_height   = TANKS.singular.fuel.h;
result.tanks.total_mass  = TANKS.singular.total_mass;
result.tanks.radius      = TANKS.singular.radius;

result.copv.mass          = IN.press.copv.mass;
result.copv.gas_mass      = IN.press.copv.gas_mass;
result.copv.amount        = IN.press.copv.amount;
result.copv.volume        = IN.press.copv.max_volume;
result.copv.max_press_psi = IN.press.copv.max_pressure / 6894.76;

result.structures.mass = OUT.Vehicle.MassDistribution.Structures;

result.avionics.battery_capacity = AVI.battery_capacity;
result.avionics.num_cells        = AVI.num_cells;
result.avionics.mass             = AVI.mass;

result.engine.eta_cstar   = IN.propulsion.eta_cstar;
result.engine.throat_area = IN.propulsion.At / 0.00064516;
result.engine.eps         = IN.propulsion.eps;
result.engine.inj_mass    = ENGINE.inj_mass;
result.engine.TCA_mass    = ENGINE.TCA_mass;
result.engine.TVC_mass    = ENGINE.TVC_mass;
result.engine.total_mass  = ENGINE.mass;

% altitude    = sim_out.altitude.Data;
% velocity    = sim_out.velocity.Data;
% 
% result.flight.max_altitude    = max(altitude);
% result.flight.max_ascent_vel  = max(velocity);
% result.flight.max_descent_vel = abs(min(velocity));
% result.flight.landing_vel     = abs(velocity(end));
% result.flight.time            = sim_out.altitude.Time(end);

x_pos = sim_out.x.Data;
y_pos = sim_out.y.Data;
z_pos = sim_out.z.Data; %add negative during plotting not now
altitude   = -(sim_out.z.Data); %plot this directly
ox_mass  = sim_out.ox_mass.Data;
fu_mass  = sim_out.fuel_mass.Data;
cg       = sim_out.cg.Data;
t        = sim_out.x.Time;
% Vertical velocity (derivative of z)
vz       = gradient(altitude, t);
thrust_data = sim_out.thrust.Data;
vel_data = sim_out.velocity.Data;
pitch_data = sim_out.pitch.Data;
roll_data = sim_out.roll.Data;
yaw_data = sim_out.yaw.Data;


result.flight.max_altitude    = max(altitude);
result.flight.max_ascent_vel  = max(vz);
result.flight.max_descent_vel = abs(min(vz));
result.flight.landing_vel     = abs(vz(end));
result.flight.flight_time     = t(end);
result.flight.final_ox_mass   = ox_mass(end);
result.flight.final_fu_mass   = fu_mass(end);
result.flight.cg_init         = cg(1);
result.flight.cg_final        = cg(end);

% Full time series for plotting
result.timeseries.t        = t;
result.timeseries.altitude = altitude;
result.timeseries.x        = x_pos;
result.timeseries.y        = y_pos;
result.timeseries.z        = z_pos;
result.timeseries.ox_mass  = ox_mass;
result.timeseries.fu_mass  = fu_mass;
result.timeseries.tot_mass = fu_mass + ox_mass;
result.timeseries.vz       = vz;
result.timeseries.thrust   = thrust_data;
result.timeseries.vx       = vel_data(:,1);
result.timeseries.vy       = vel_data(:,2);
result.timeseries.vz_raw   = vel_data(:,3);
result.timeseries.roll     = roll_data;
result.timeseries.pitch    = pitch_data;
result.timeseries.yaw      = yaw_data;



result.status.success = true;
result.status.error   = '';


% Sanity check — flag diverged runs
if any(isinf(altitude)) || any(isnan(altitude)) || max(altitude) > 500
    result.status.success = false;
    result.status.error   = 'Simulation diverged - states went to Inf/NaN';
    return
end

end
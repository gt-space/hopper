function result = mc_main_ctv(scenario, sim_out)
%MC_MAIN_CTV Collapse one CTV Monte Carlo run into a result struct.
%   Same result struct as sims/hopper_6dof_sizer/mc_main.m, so mc_plots_ctv
%   works on it unchanged. Unlike mc_main it does not re-run the sizing:
%   mc_sim_setup_ctv has already put IN/VEH/TANKS/STRUCT/OUT/AVI/ENGINE in the
%   base workspace, and re-running CTV's sizing here would just recompute them.

for v = {'IN', 'VEH', 'TANKS', 'STRUCT', 'OUT', 'AVI', 'ENGINE'}
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', v{1}))
        error('mc_main_ctv:missingBaseVars', ...
            'Base workspace has no %s - run mc_sim_setup_ctv first.', v{1});
    end
end
IN     = evalin('base', 'IN');
VEH    = evalin('base', 'VEH');     %#ok<NASGU>
TANKS  = evalin('base', 'TANKS');
STRUCT = evalin('base', 'STRUCT');  %#ok<NASGU>
OUT    = evalin('base', 'OUT');
AVI    = evalin('base', 'AVI');
ENGINE = evalin('base', 'ENGINE');

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
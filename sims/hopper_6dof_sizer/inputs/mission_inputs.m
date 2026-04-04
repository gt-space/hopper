function IN = mission_inputs(scenario_params)
if nargin < 1; scenario_params = []; end

%  GLOBAL CONSTANTS
IN.const.g0 = 9.80665;              % m/s^2
IN.const.R_univ = 8.314462618;      % J/mol-K
IN.const.FS_struct = 2;           % structural factor of safety

%  MISSION CONSTRAINTS (checked after flight simulation)
IN.mission.target_altitude = 51;    % m
IN.mission.hover_time_55m = 2;      % s
IN.mission.max_ascent_vel = 7;     % m/s
IN.mission.max_descent_vel = 7;    % m/s
IN.mission.hover_time_1m = 1;       %s
IN.mission.landing_vel = 0;         % m/s at 1 m
IN.mission.landing_alt = 1;         % m
IN.mission.TWR_hover_landing = 1;
IN.mission.max_vehicle_height = 2.5; % m

%  MASS GROWTH MARGIN
IN.margins.mass_growth = 20 * 0.453592;  % lb -> kg

%  STRUCTURES
IN.structures.payload_mass = 15; % kg, CPLC Requirement
IN.structures.extra_payload_mass = 16.13; % kg
IN.structures.payload_cg_z = [];            % m (from cad/aidan)

%  AVIONICS
IN.avionics.standby_hours = 2.0; % hr
IN.avionics.flight_time = 2 * 30; % sec
IN.avionics.lcoms_time = 10 * 60; % min --> sec
IN.avionics.battery_voltage = 48; % V
IN.avionics.boards.voltage = 12; % V
IN.avionics.boards.current = 0.5; % A

%  PROPULSION
IN.propulsion.oxidizer = 'Oxygen';
IN.propulsion.fuel = 'Dodecane'; %AKA KEROSENE
IN.propulsion.oxidizer_mass = 15; % kg -> internal
IN.propulsion.fuel_mass = 10; % kg -> internal
IN.propulsion.ox_vol = 0.020; % m^3
IN.propulsion.fu_vol = 0.0152; % m^3
IN.propulsion.ox_press = 500 * 6894.76; % Pa
IN.propulsion.fu_press = 500 * 6894.76; % Pa
IN.propulsion.fu_temp = 293; % K
IN.propulsion.ox_temp = 93; % K
IN.propulsion.ox_sys_CdA = 1.7827E-5; % m^2
IN.propulsion.fu_sys_CdA = 1.1E-5; % m^2
IN.propulsion.eta_cstar = 0.85;
IN.propulsion.At = 1.138 * 0.00064516; % in^2 --> m^2
IN.propulsion.eps = 3.655;
IN.propulsion.nominal_thrust = 500 * 4.44822; % N
IN.propulsion.throttle_range = [0.4 1.1];
IN.solenoid_valve_quantity = 5;
IN.valves.voltage = 12; % V
IN.valves.current = 1; % A

% ENGINE
IN.inj_material = "SS316L";
IN.TCA_material = "AlSi10Mg";
IN.off_axis_y = 0; % m^2  % Body Frame
IN.off_axis_x = 0; % m^2 

%  PRESSURIZATION MODE
IN.press.mode = 'copv';  % 'copv' or 'autogenous'

% --- COPV Parameters ---
% subscale COPV params
IN.press.copv.max_pressure = 4500 * 6894.76; % Pa
IN.press.copv.max_volume = 0.0267;               % m^3
IN.press.copv.gas = 'GN2';
IN.press.copv.initial_temperature = 293;             % K
IN.press.copv.mass = 5; % kg
IN.press.copv.amount = 3;
IN.press.copv.gas_mass = 8.27;

% --- Autogenous Parameters ---
IN.press.autogenous.gas = 'N2O';
IN.press.autogenous.initial_temperature = 300;       % K
IN.press.autogenous.piston_friction = 0;      % Pa

%  TANK OPTIONS
IN.tanks.geometry = 'stacked'; % stacked, concentric, clustered
IN.tanks.material = 'Al6061';
IN.tanks.max_radius = 4*0.0254; % in -> m
IN.tanks.lateral_damping = 0.02;
IN.tanks.axial_damping = 0.1;

%  TVC
IN.tvc.max_gimbal_angle = deg2rad(15); % ?
IN.tvc.max_gimbal_rate = deg2rad(30); % rad/s ?
IN.tvc.voltage = 24; % V
IN.tvc.current = 3; % A
IN.tvc.actuator_latency = 0.01; % s
IN.tvc.resolution = deg2rad(0.1);
IN.tvc.damping_ratio = 0.5; 

% RCS
IN.rcs.max_thrust = 10; % N
IN.rcs.throttle_rate = 25; % N/s
IN.rcs.max_speed = 62500; % RPM
IN.rcs.min_speed = 0.05; % RPM
IN.rcs_time_constant = 0.1;
IN.rcs.resolution = 2.5;

% THROTTLE VALVE
IN.throttle_valve.voltage = 24; % V
IN.throttle_valve.current = 1; % A
IN.throttle_valve.mass = 2.3; % kg
IN.throttle_valve.rate_limit = 90 * 4.44822; % N
IN.throttle_valve_latency = 0.01; % s
IN.throttle_valve.resolution = 5; % N

%  LANDING LEGS
IN.legs.material = 'Al6061';
IN.legs.tip_factor = 2; % accounts for non-symmetric landing (leg takes 50% extra load)
IN.legs.feet_radial_dist = 0.5; % m

%cg values
IN.engine_cg = 0.3556;
IN.mount_cg = 0.45;
IN.ox_tank_cg = 1.854;
IN.ox_tank_wall_thick = 0.003175;
IN.fu_tank_cg = 1.0668;
IN.fu_tank_wall_thick = 0.003175;
IN.structures_cg = 1.193;
IN.avi_cg = 1.4224;
IN.fluids_cg = 1.0668;
IN.payload_cg = 2;
IN.copv_r = 0.18;
IN.copv_t = 13.21/1000;
IN.copv_h = 0.56;
IN.copv_zcg = 1.21;
IN.copv_ycg = 0.23;


%  MODAL PARAMETERS
IN.vehicle_damping_ratio = 0.05; 
IN.vehicle_natural_frequency = 10; % Hz
%=======
% Wind Defaults
% 1 - Monte Carlo , 2 - constant , 3 - 4D Array Historical Wind
IN.wind.mode = 2; %
IN.wind.uwind = 2.5; %default values for uwnd/vwnd with NO Monte Carlo
IN.wind.vwind = 2.5;
%>>>>>>> Stashed changes

%  SIMULATION SETTINGS
IN.sim.dt = 0.01;                  % s
IN.sim.t_max = 60;                 % s
IN.sim.use_simulink = true;


% Monte Carlo scenario struct creation
if  ~isempty(scenario_params)
    %IN.propulsion.oxidizer_mass       = scenario_params.ox_mass;
    IN.propulsion.fuel_mass           = scenario_params.fuel_mass;
    IN.propulsion.eta_cstar           = scenario_params.cstar;
    IN.mass_factor                    = scenario_params.mass_factor;
    IN.tanks.lateral_damping          = scenario_params.slosh_lateral_damping;
    IN.tanks.axial_damping            = scenario_params.slosh_axial_damping;
    IN.throttle_valve.rate_limit      = scenario_params.throttle_rate_limit * 4.44822;
    IN.throttle_valve_latency         = scenario_params.throttle_latency;
    IN.tvc.max_gimbal_rate            = deg2rad(scenario_params.tvc_actuator_rate_limit);
    IN.tvc.actuator_latency           = scenario_params.tvc_actuator_latency;
    IN.off_axis_y                     = scenario_params.engine_off_axis_y;
    IN.off_axis_x                     = scenario_params.engine_off_axis_z;
    IN.wind.uwind                     = scenario_params.uwind;
    IN.wind.vwind                     = scenario_params.vwind;
    IN.wind.mode                      = 1; 
end

end


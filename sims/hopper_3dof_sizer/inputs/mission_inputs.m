function IN = mission_inputs()
%  GLOBAL CONSTANTS
IN.const.g0 = 9.80665;              % m/s^2
IN.const.R_univ = 8.314462618;      % J/mol-K
IN.const.FS_struct = 1.5;           % structural factor of safety

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
IN.structures.extra_payload_mass = 44.13; % kg
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
IN.propulsion.fuel = 'Kerosene';
IN.propulsion.oxidizer_mass = 19.43; % kg -> internal
IN.propulsion.fuel_mass = 9.7; % kg -> internal
IN.propulsion.ox_vol = 0.020; % m^3
IN.propulsion.fu_vol = 0.0152; % m^3
IN.propulsion.nominal_thrust = 500 * 4.44822; % N
IN.propulsion.throttle_range = [0.4 1.1];
IN.solenoid_valve_quantity = 5;
IN.valves.voltage = 12; % V
IN.valves.current = 1; % A

% ENGINE
IN.inj_material = "SS316L";
IN.TCA_material = "AlSi10Mg";

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

%  TVC
IN.tvc.max_gimbal_angle = deg2rad(5); % ?
IN.tvc.max_gimbal_rate = deg2rad(20); % rad/s ?
IN.tvc.voltage = 24; % V
IN.tvc.current = 3; % A

% THROTTLE VALVE
IN.throttle_valve.voltage = 24; % V
IN.throttle_valve.current = 1; % A
IN.throttle_valve.mass = 2.3; % kg

%  LANDING LEGS
IN.legs.material = 'Al6061';
IN.legs.tip_factor = 1.5; % accounts for non-symmetric landing (leg takes 50% extra load)

%  SIMULATION SETTINGS
IN.sim.dt = 0.01;                  % s
IN.sim.t_max = 60;                 % s
IN.sim.use_simulink = true;

end
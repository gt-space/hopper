function IN = mission_inputs(scenario_params)
if nargin < 1; scenario_params = []; end

%  GLOBAL CONSTANTS
IN.const.g0 = 9.80665;              % m/s^2
IN.const.R_univ = 8.314462618;      % J/mol-K
IN.const.FS_struct = 2;             % structural factor of safety

%  MISSION CONSTRAINTS (checked after flight simulation)
% checkConstraints judges every flight against the fields marked [req];
% change a requirement here and it applies everywhere (MC summaries,
% meets_requirements). CI's looser regression limits: tests/montecarlo/mcLimits.m
IN.mission.target_altitude = 51;    % m          [req] minimum peak altitude
IN.mission.hover_time_55m = 2;      % s
IN.mission.three_sec_check = false; % true: hold 3 s on the pad before ascent (ThreeSecTrajectory2)
IN.mission.max_ascent_vel = 7;      % m/s        [req]
IN.mission.max_descent_vel = 7;     % m/s        [req]
IN.mission.max_impact_vel = 0.5;    % m/s        [req] vertical speed at ground contact
IN.mission.max_landing_offset = 0.5; % m         [req] from pad center (pad radius)
IN.mission.max_tilt_deg = 10;       % deg        [req] body axis from vertical, whole flight
IN.mission.min_twr = 1.0;           % -          [req] initial TWR at full throttle
IN.mission.hover_time_1m = 1;       % s
IN.mission.landing_vel = 0;         % m/s at 1 m
IN.mission.landing_alt = 1;         % m
IN.mission.TWR_hover_landing = 1;
IN.mission.max_vehicle_height = 2.44; % m

%  STRUCTURES
IN.structures.payload_mass = 15;    % kg, CPLC Requirement
IN.structures.structures_mass = 34.26; % kg from CAD

%  AVIONICS
IN.avionics.standby_hours = 2.0;    % hr
IN.avionics.flight_time = 2 * 30;   % sec
IN.avionics.lcoms_time = 10 * 60;   % min --> sec
IN.avionics.battery_voltage = 48;   % V
IN.avionics.boards.voltage = 12;    % V
IN.avionics.boards.current = 0.5;   % A
IN.avionics.avi_mass = 11.52;       % kg from CAD

%  PROPULSION
IN.propulsion.oxidizer = 'O2';      % Oxygen - coolprop
IN.propulsion.fuel = 'Dodecane';    % AKA KEROSENE - coolprop
IN.propulsion.oxidizer_mass = 16;   % kg -> internal
IN.propulsion.fuel_mass = 13;       % kg -> internal
IN.propulsion.ox_vol = 0.01638;     % m^3
IN.propulsion.fu_vol = 0.02034;     % m^3
IN.propulsion.ox_press = 500 * 6894.76; % Pa
IN.propulsion.fu_press = 500 * 6894.76; % Pa
IN.propulsion.fu_temp = 293;        % K
IN.propulsion.ox_temp = 93;         % K
IN.propulsion.ox_sys_CdA = 1.7827E-5; % m^2
IN.propulsion.fu_sys_CdA = 1.1E-5;  % m^2
IN.propulsion.eta_cstar = 0.85;
IN.propulsion.At = 1.056 * 0.00064516; % in^2 --> m^2
IN.propulsion.eps = 2.93;
IN.propulsion.nominal_thrust = 500 * 4.44822; % N
IN.propulsion.throttle_range = [0.4 1.1];
IN.solenoid_valve_quantity = 5;
IN.valves.voltage = 12;             % V
IN.valves.current = 1;              % A
IN.fluids_mass = 28.41;             % kg from CAD

% ENGINE
IN.inj_material = "SS316L";
IN.TCA_material = "AlSi10Mg";
IN.off_axis_y = 0;                  % m^2  % Body Frame
IN.off_axis_x = 0;                  % m^2 
IN.engine_mass = 37.83;             % kg from CAD

%  PRESSURIZATION MODE
IN.press.mode = 'copv';             % 'copv' or 'autogenous'

% --- COPV Parameters ---
IN.press.copv.max_pressure = 4500 * 6894.76; % Pa
IN.press.copv.max_volume = 0.0267;  % m^3
IN.press.copv.gas = 'GN2';
IN.press.copv.initial_temperature = 293; % K
IN.press.copv.mass = 11;            % kg
IN.press.copv.amount = 3;
IN.press.copv.gas_mass = 6;         % kg
IN.press.copv.r = 0.18;             % m
IN.press.copv.t = 13.21/1000;       % m
IN.press.copv.h = 0.56;             % m

% --- Autogenous Parameters ---
IN.press.autogenous.gas = 'N2O';
IN.press.autogenous.initial_temperature = 300; % K
IN.press.autogenous.piston_friction = 0; % Pa

%  TANK OPTIONS
IN.tanks.geometry = 'stacked';      % stacked, concentric, clustered
IN.tanks.material = 'Al6061';
IN.tanks.max_radius = 4.25*0.0254;  % in -> m
% As-built tank lengths. These override the sizer TANKS.singular.*.h for
% CG/MOI and slosh. The Slosh and CG // MOI
% r_*_tank constants in hopper_6dof_NED_v2.slx use the same values.
IN.tanks.ox_h_as_built = 29.18*0.0254;  % LOx, in -> m
IN.tanks.fu_h_as_built = 18.91*0.0254;  % RP-1, in -> m
IN.tanks.lateral_damping = 0.08;   % <-- MUST BE UNCOMMENTED
IN.tanks.axial_damping = 0.3;      % <-- MUST BE UNCOMMENTED
IN.tanks.damping_factor_ox = 1.0;  % MC: scale on the baffle damping table (LOx)
IN.tanks.damping_factor_fu = 1.0;  % MC: scale on the baffle damping table (RP-1)

%  TVC
IN.tvc.max_gimbal_angle = deg2rad(15); % deg
IN.tvc.max_gimbal_rate = deg2rad(15);  % deg
IN.tvc.voltage = 24;                % V
IN.tvc.current = 3;                 % A
IN.tvc.actuator_latency = 0.01;     % s
IN.tvc.rise_time = 0.05;            % s
IN.tvc.resolution = deg2rad(0.1);   % deg
IN.tvc.damping_ratio = 0.5; 

% RCS
IN.rcs.max_thrust = 16;             % N
IN.rcs.throttle_rate = 50;          % N/s
IN.rcs_time_constant = 0.1;
IN.rcs.resolution = 0.1;            % N

% THROTTLE VALVE
IN.throttle_valve.voltage = 24;     % V
IN.throttle_valve.current = 1;      % A
IN.throttle_valve.mass = 2.3;       % kg
IN.throttle_valve.rate_limit = 180 * 4.44822; % N
IN.throttle_valve_latency = 0.01;   % s
IN.throttle_valve.resolution = 5;   % N

%  LANDING LEGS
IN.legs.material = 'Al6061';
IN.legs.tip_factor = 2;             % non-symmetric landing factor
IN.legs.feet_radial_dist = 0.5;     % m
IN.legs.damping_const = 1.0;        % strut damping

% CG VALUES
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

% MODAL PARAMETERS
IN.vehicle_damping_ratio = 0.01; 
IN.vehicle_natural_frequency = 10;  % Hz

% WIND DEFAULTS
IN.wind.mode = 2; 
IN.wind.uwind = 10; 
IN.wind.vwind = 10;

% MASS FACTOR
IN.mass_factor = 1.0;               % default / nominal case
IN.cg_factor = 1.0;                 % dry CG axial position factor (MC)
IN.inertia_factor = 1.0;            % dry inertia factor (MC)

% MONTE CARLO SCENARIO STRUCT CREATION (SAFEGUARDED)
if ~isempty(scenario_params) && isfield(scenario_params, 'ox_mass')
    IN.propulsion.oxidizer_mass       = scenario_params.ox_mass;
    IN.propulsion.fuel_mass           = scenario_params.fuel_mass;
    IN.propulsion.eta_cstar           = scenario_params.cstar;
    IN.mass_factor                    = scenario_params.mass_factor;
    if isfield(scenario_params, 'slosh_lateral_damping'), IN.tanks.lateral_damping = scenario_params.slosh_lateral_damping; end
    if isfield(scenario_params, 'slosh_axial_damping'),   IN.tanks.axial_damping   = scenario_params.slosh_axial_damping;   end
    % slosh damping dispersion: log-uniform scale on the baffle damping table
    if isfield(scenario_params, 'log2_slosh_damping_ox'), IN.tanks.damping_factor_ox = 2^scenario_params.log2_slosh_damping_ox; end
    if isfield(scenario_params, 'log2_slosh_damping_fu'), IN.tanks.damping_factor_fu = 2^scenario_params.log2_slosh_damping_fu; end
    IN.throttle_valve.rate_limit      = scenario_params.throttle_rate_limit * 4.44822;
    IN.throttle_valve_latency         = scenario_params.throttle_latency;
    IN.tvc.max_gimbal_rate            = deg2rad(scenario_params.tvc_actuator_rate_limit);
    IN.tvc.actuator_latency           = scenario_params.tvc_actuator_latency;
    IN.off_axis_y                     = scenario_params.engine_off_axis_y;
    IN.off_axis_x                     = scenario_params.engine_off_axis_z;
    IN.wind.uwind                     = scenario_params.uwind;
    IN.wind.vwind                     = scenario_params.vwind;
    IN.wind.mode                      = 1; 
    if isfield(scenario_params, 'cg_factor'),      IN.cg_factor      = scenario_params.cg_factor;      end
    if isfield(scenario_params, 'inertia_factor'), IN.inertia_factor = scenario_params.inertia_factor; end
end
end
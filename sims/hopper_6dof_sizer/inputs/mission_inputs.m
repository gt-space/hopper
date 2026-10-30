function IN = mission_inputs(scenario_params)
% Mission inputs: requirements, trajectory, environment, sizer settings and
% Monte Carlo dispersions. The vehicle itself (geometry, masses, hardware)
% comes from its config file, vehicles/<vehicle>/<vehicle>_config.m, for
% the vehicle hopperVehicle selects. Another vehicle's mission overrides
% are in inputs/<vehicle>_inputs.m, applied last.
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
% Landing limits follow the GNC PDR (slide 107): radius < 1 m, normal
% velocity < 1.5 m/s, lateral velocity < 0.5 m/s, touchdown tilt < 2 deg.
IN.mission.max_ascent_vel = 10;     % m/s        [req]
IN.mission.max_descent_vel = 10;    % m/s        [req]
IN.mission.max_impact_vel = 1.5;    % m/s        [req] vertical (normal) speed at ground contact
IN.mission.max_lateral_vel = 0.5;   % m/s        [req] horizontal speed at ground contact
IN.mission.max_landing_offset = 1;  % m          [req] from pad center
IN.mission.max_touchdown_tilt_deg = 2; % deg     [req] body axis from vertical at ground contact
IN.mission.max_tilt_deg = 10;       % deg        [req] body axis from vertical, whole flight
IN.mission.min_twr = 1.0;           % -          [req] initial TWR at full throttle
IN.mission.hover_time_1m = 1;       % s
IN.mission.landing_vel = 0;         % m/s at 1 m
IN.mission.landing_alt = 1;         % m
IN.mission.engine_cutoff_alt = 0.05; % m, engine off below this height while descending (Landing_Engine_Cutoff)
IN.mission.TWR_hover_landing = 1;
% Reference trajectory (Trajectory2): minimum-jerk climb to apex, hover,
% minimum-jerk descent to z_td, then a constant-sink segment to touchdown
IN.mission.trajectory.apex = 52;    % m
IN.mission.trajectory.t_up = 12;    % s
IN.mission.trajectory.t_hover = 2;  % s
IN.mission.trajectory.t_down = 13;  % s
IN.mission.trajectory.t_term = 2.5; % s, constant-sink terminal segment
IN.mission.trajectory.z_td = 0.4;   % m, where the constant-sink segment starts
IN.mission.trajectory.v_td = 0.3;   % m/s, touchdown sink rate
IN.mission.trajectory.burn_mass = 21; % kg, propellant the thrust reference assumes is burned

%  AVIONICS
IN.avionics.standby_hours = 2.0;    % hr
IN.avionics.flight_time = 2 * 30;   % sec
IN.avionics.lcoms_time = 10 * 60;   % min --> sec
IN.avionics.battery_voltage = 48;   % V
IN.avionics.boards.voltage = 12;    % V
IN.avionics.boards.current = 0.5;   % A

%  PROPULSION
IN.propulsion.oxidizer = 'O2';      % Oxygen - coolprop
IN.propulsion.fuel = 'Dodecane';    % AKA KEROSENE - coolprop
IN.propulsion.ox_press = 500 * 6894.76; % Pa
IN.propulsion.fu_press = 500 * 6894.76; % Pa
IN.propulsion.fu_temp = 293;        % K
IN.propulsion.ox_temp = 93;         % K
IN.propulsion.ox_sys_CdA = 1.7827E-5; % m^2
IN.propulsion.fu_sys_CdA = 1.1E-5;  % m^2
IN.propulsion.eta_cstar = 0.85;
IN.propulsion.At = 1.056 * 0.00064516; % in^2 --> m^2
IN.propulsion.eps = 2.93;
IN.solenoid_valve_quantity = 5;
IN.valves.voltage = 12;             % V
IN.valves.current = 1;              % A

% ENGINE
IN.inj_material = "SS316L";
IN.TCA_material = "AlSi10Mg";

%  PRESSURIZATION MODE
IN.press.mode = 'copv';             % 'copv' or 'autogenous'

% --- COPV Parameters ---
IN.press.copv.max_pressure = 4500 * 6894.76; % Pa
IN.press.copv.max_volume = 0.0267;  % m^3
IN.press.copv.gas = 'GN2';
IN.press.copv.initial_temperature = 293; % K

% --- Autogenous Parameters ---
IN.press.autogenous.gas = 'N2O';
IN.press.autogenous.initial_temperature = 300; % K
IN.press.autogenous.piston_friction = 0; % Pa

%  TANK OPTIONS
IN.tanks.lateral_damping = 0.08;   % <-- MUST BE UNCOMMENTED
IN.tanks.axial_damping = 0.3;      % <-- MUST BE UNCOMMENTED
IN.tanks.damping_factor_ox = 1.0;  % MC: scale on the baffle damping table (LOx)
IN.tanks.damping_factor_fu = 1.0;  % MC: scale on the baffle damping table (RP-1)

%  TVC
IN.tvc.voltage = 24;                % V
IN.tvc.current = 3;                 % A

% THROTTLE VALVE
IN.throttle_valve.voltage = 24;     % V
IN.throttle_valve.current = 1;      % A
IN.throttle_valve.mass = 2.3;       % kg

%  LANDING LEGS
IN.legs.material = 'Al6061';
IN.legs.tip_factor = 2;             % non-symmetric landing factor

% WIND DEFAULTS
IN.wind.mode = 2;
IN.wind.uwind = 10;
IN.wind.vwind = 10;

% AERODYNAMICS (Wind_Forces)
IN.aero.rho = 1.223;                % kg/m^3

% VEHICLE: geometry, masses and hardware come from the vehicle's config
% file (vehicles/<vehicle>/<vehicle>_config.m, docs/VEHICLE_CONFIG.md)
IN = vehicleInputs(IN, vehicleConfig(hopperVehicle()));

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

% Other vehicles: their overrides win over the hopper values above
switch hopperVehicle()
    case 'ctv', IN = ctv_inputs(IN, scenario_params);
end
end
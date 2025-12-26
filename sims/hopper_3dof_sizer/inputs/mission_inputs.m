function IN = mission_inputs()

%% =======================
%  GLOBAL CONSTANTS
%% =======================
IN.const.g0 = 9.80665;              % m/s^2
IN.const.R_univ = 8.314462618;      % J/mol-K
IN.const.FS_struct = 1.5;           % structural factor of safety

%% =======================
%  MISSION REQUIREMENTS
%% =======================
IN.mission.target_altitude = 55;    % m
IN.mission.hover_time = 5;          % s
IN.mission.max_ascent_vel = 20;     % m/s
IN.mission.landing_vel = 0;         % m/s at 1 m
IN.mission.landing_alt = 1;         % m
IN.mission.gravity_loss_factor = 1.05;

%% =======================
%  MASS GROWTH MARGIN
%% =======================
IN.margins.mass_growth = 20 * 0.453592;  % kg

%% =======================
%  STRUCTURES
%% =======================
IN.structures.payload_mass = 34 * 0.453592; % kg
IN.structures.payload_cg_z = [];            % m (filled later)

%% =======================
%  AVIONICS
%% =======================
IN.avionics.mass = [];              % kg
IN.avionics.power = [];             % W
IN.avionics.hours = 1.0;            % hr

%% =======================
%  PROPULSION
%% =======================
IN.propellant.oxidizer = 'N2O';
IN.propellant.fuel = 'IPA';

IN.propulsion.nominal_thrust = 500 * 4.44822; % N
IN.propulsion.throttle_range = [0.3 1.0];
IN.propulsion.Isp_vac = [];         % s
IN.propulsion.OF = [];              % mixture ratio
IN.propulsion.chamber_pressure = []; % Pa
IN.propulsion.injector_dp_frac = 0.2; % ΔP/Pc

%% =======================
%  PRESSURIZATION MODE
%% =======================
IN.press.mode = 'copv';  % 'copv' or 'autogenous'

% --- COPV Parameters ---
IN.press.copv.max_pressure = 4500 * 6894.76; % Pa
IN.press.copv.max_volume = [];               % m^3
IN.press.copv.gas = 'GN2';
IN.press.copv.temperature = 300;             % K

% --- Autogenous Parameters ---
IN.press.autogenous.gas = 'N2O';
IN.press.autogenous.temperature = 300;       % K
IN.press.autogenous.piston_friction = 0;      % Pa

%% =======================
%  TANK OPTIONS
%% =======================
IN.tanks.geometry = 'cylindrical'; % future-proof
IN.tanks.material = 'Al2219';
IN.tanks.allow_clustered = true;
IN.tanks.allow_concentric = true;
IN.tanks.max_radius = [];          % m

%% =======================
%  TVC
%% =======================
IN.tvc.max_gimbal_angle = deg2rad(5);
IN.tvc.max_gimbal_rate = deg2rad(20); % rad/s

%% =======================
%  LANDING LEGS
%% =======================
IN.legs.material = 'Al6061';
IN.legs.tip_factor = 1.5;

%% =======================
%  BATTERIES
%% =======================
IN.battery.energy_density = 200;   % Wh/kg
IN.battery.margin = 1.5;

%% =======================
%  SIMULATION SETTINGS
%% =======================
IN.sim.dt = 0.01;                  % s
IN.sim.t_max = 60;                 % s
IN.sim.use_simulink = true;

end
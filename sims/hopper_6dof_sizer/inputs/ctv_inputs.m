function IN = ctv_inputs(IN, scenario)
%CTV_INPUTS  Controls test vehicle: overrides on top of the hopper inputs.
%   IN = ctv_inputs(IN)            what mission_inputs returns when
%                                  hopperVehicle is 'ctv'
%                                  (variantConfig('VEHICLE', 'ctv'))
%   IN = ctv_inputs(IN, scenario)  one Monte Carlo scenario
%                                  (inputs/mc_params_ctv.json), applied
%                                  after the CTV values
%
%   The CTV is a ~2.6 kg electric vehicle: one gimbaled EDF (~35 N) for
%   thrust and pitch/yaw, four small props for roll. It flies the same
%   actuator model as the hopper and burns no propellant, so it runs with
%   SLOSH=off and MASS_PROPS=constant (variantConfig's CTV defaults).
%
%   Mass, CG and inertia come from CAD (as in ctv/inputs/mission_inputs.m).
%   Values marked ASSUMPTION are placeholders until measured; check them
%   before trusting a CTV result.

% ---- CAD mass properties (CAD frame: origin at the base centre, Y up) --
IN.ctv.mass = 2614.88 / 1000;                       % kg
IN.ctv.cg_cad = [-41.76; 349.21; 6.07] * 1e-3;      % m, CAD X Y Z
IN.ctv.I_cad = [103179574.33   1374767.18     68051.49;
                  1374767.18  12321923.43  -1150369.65;
                    68051.49  -1150369.65  103816055.97] * 1e-9;  % kg m^2 about the CG, CAD axes
IN.ctv.thrust = 7.875 * 4.44822;                    % N, EDF static thrust

% Body x is the thrust axis (CAD Y); y, z follow right-handed (CAD Z, CAD X)
P = [0 1 0; 0 0 1; 1 0 0];                          % body = P * CAD
cg_body = P * IN.ctv.cg_cad;
IN.mass.dry_cg = cg_body(1);                        % m, axial CG from the base
IN.ctv.I_body = P * IN.ctv.I_cad * P';             % kg m^2, body axes, full CAD tensor
IN.ctv.cg_lateral = cg_body(2:3);                   % m; the sim uses the axial CG only
% Principal moments only, as with the CG: the CAD products of inertia
% (up to 9 % of the small roll inertia) make the LQR's flight drift 8.5 m
% and tilt 20 deg in calm air, with no disturbance at all. Investigate that
% coupling before trusting a CTV flight with the full tensor.
IN.mass.dry_I = diag(diag(IN.ctv.I_body));          % kg m^2

% ---- Propellant: none. The propellant model needs tanks that aren't
% empty (the depletion cutoff stops thrust at zero), so it carries 1 g
% placeholders inside the CAD mass; MASS_PROPS=constant keeps them fixed.
IN.propulsion.burns_propellant = false;
IN.propulsion.oxidizer_mass = 1e-3;                 % kg, placeholder
IN.propulsion.fuel_mass = 1e-3;                     % kg, placeholder
IN.mass.dry_override = IN.ctv.mass - IN.propulsion.oxidizer_mass - IN.propulsion.fuel_mass;

% ---- EDF and gimbal ----
IN.propulsion.nominal_thrust = IN.ctv.thrust;      % N
IN.propulsion.throttle_range = [0.2 1.0];          % ASSUMPTION: EDF idle to full
IN.tvc.max_gimbal_angle = deg2rad(8);              % from ctv/inputs/mission_inputs.m
IN.engine_mass = 0.35;                             % kg, ASSUMPTION: EDF + gimbal (inside the CAD mass)
IN.engine_cg = 0.08;                               % m, ASSUMPTION: EDF CG above the base
IN.mount_cg = 0.15;                                % m, ASSUMPTION: gimbal pivot above the base
IN.mission.max_vehicle_height = 0.75;              % m, ASSUMPTION

% EDF throttle: an ESC, not the hopper's throttle valve (5 N steps there
% are 20 % of the CTV's weight and made thrust bang-bang)
IN.throttle_valve.resolution = 0.05;               % N, ASSUMPTION
IN.tvc.max_gimbal_rate = deg2rad(200);             % rad/s, ASSUMPTION: hobby servo (hopper 15 deg/s)

% ---- Roll props ----
% Four 65 mm props on coreless motors (ctv/inputs/mission_inputs.m), still
% at that file's placeholder thrust. TODO: datasheet value. The roll
% authority matters: with the current LQR the CTV holds attitude under a
% 1 mm engine offset or 1.5 m/s wind at 16 N but tumbles at 2 N or less
% (typical 65 mm props give ~0.4 N), so check this before flying.
IN.rcs.max_thrust = 16;                            % N, placeholder (ctv/inputs/mission_inputs.m)

% ---- LQR weights: the hopper's (LinerizationMaster). Set CTV weights with
% IN.gnc.lqr.Q_diag / R_diag / Sf_diag (13 / 4 / 13 diagonal entries).
% Tried 2026-10-28: lateral position x10 and velocity x30 landed one windy
% case within 0.25 m but did worse than the hopper's weights across a
% 10-run Monte Carlo, so the CTV still needs real tuning (docs/VARIANTS.md).

% ---- Landing legs ----
IN.legs.feet_radial_dist = 0.15;                   % m, ASSUMPTION

% ---- Aerodynamics (Wind_Forces) ----
IN.aero.area = 0.15;                               % m^2, ASSUMPTION: side area
IN.aero.cd = 1.2;                                  % ASSUMPTION
IN.aero.cp = 0.35;                                 % m, ASSUMPTION: centre of pressure from the base

% ---- Wind: calm. ASSUMPTION: the CTV flies in calm air; the hopper's
% 10 m/s default puts ~11 N of side force (~40 % of weight) on it
IN.wind.uwind = 0;                                 % m/s
IN.wind.vwind = 0;                                 % m/s

% ---- Mission: 4 m hop with limits scaled to the CTV ----
IN.mission.target_altitude = 4;                    % m  [req]
IN.mission.max_ascent_vel = 2;                     % m/s [req] ASSUMPTION (CTV inputs had 26)
IN.mission.max_descent_vel = 2;                    % m/s [req] ASSUMPTION
IN.mission.trajectory.apex = 4.5;                  % m
IN.mission.trajectory.t_up = 6;                    % s (peak 1.4 m/s)
IN.mission.trajectory.t_hover = 2;                 % s
IN.mission.trajectory.t_down = 6;                  % s
IN.mission.trajectory.t_term = 2;                  % s
IN.mission.trajectory.z_td = 0.2;                  % m
IN.mission.trajectory.v_td = 0.1;                  % m/s
IN.mission.trajectory.burn_mass = 0;               % kg
% An EDF throttles down to the ground: cut on contact, not 5 cm above it
% (a 5 cm drop alone lands at ~1 m/s)
IN.mission.engine_cutoff_alt = 0.005;              % m

% ---- Monte Carlo dispersions (inputs/mc_params_ctv.json) ----
if nargin > 1 && ~isempty(scenario)
    s = scenario;
    if isfield(s, 'mass_factor')
        IN.mass_factor = s.mass_factor;
        IN.mass.dry_override = s.mass_factor * IN.mass.dry_override;
    end
    if isfield(s, 'cg_factor'),      IN.cg_factor = s.cg_factor; end
    if isfield(s, 'inertia_factor'), IN.inertia_factor = s.inertia_factor; end
    if isfield(s, 'thrust_factor'),  IN.propulsion.nominal_thrust = s.thrust_factor * IN.propulsion.nominal_thrust; end
    if isfield(s, 'tvc_actuator_rate_limit'), IN.tvc.max_gimbal_rate = deg2rad(s.tvc_actuator_rate_limit); end
    if isfield(s, 'tvc_actuator_latency'),    IN.tvc.actuator_latency = s.tvc_actuator_latency; end
    if isfield(s, 'engine_off_axis_y'), IN.off_axis_y = s.engine_off_axis_y; end
    if isfield(s, 'engine_off_axis_z'), IN.off_axis_x = s.engine_off_axis_z; end
    if isfield(s, 'uwind'), IN.wind.uwind = s.uwind; IN.wind.mode = 1; end
    if isfield(s, 'vwind'), IN.wind.vwind = s.vwind; IN.wind.mode = 1; end
end
end

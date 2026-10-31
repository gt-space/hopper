function IN = ctv_inputs(IN, scenario)
%CTV_INPUTS  Controls test vehicle: mission overrides on top of the hopper's.
%   IN = ctv_inputs(IN)            what mission_inputs returns when
%                                  hopperVehicle is 'ctv'
%                                  (variantConfig('VEHICLE', 'ctv'))
%   IN = ctv_inputs(IN, scenario)  one Monte Carlo scenario
%                                  (inputs/mc_params_ctv.json), applied
%                                  after the CTV values
%
%   The vehicle itself (CAD mass properties, EDF, gimbal, roll props, legs,
%   aero) is vehicles/ctv/ctv_config.m. This file has what the CTV flies:
%   the mission, wind and its Monte Carlo dispersions. Values marked
%   ASSUMPTION are placeholders until measured.

% ---- LQR weights: the hopper's (LinerizationMaster), except roll. Set CTV
% weights with IN.gnc.lqr.Q_diag / R_diag / Sf_diag (13 / 4 / 13 diagonal
% entries). Tried 2026-10-28: lateral position x10 and velocity x30 landed
% one windy case within 0.25 m but did worse than the hopper's weights
% across a 10-run Monte Carlo, so the CTV still needs real tuning
% (docs/VARIANTS.md).
% Roll: the LQR's roll-rate gain is sqrt(Q_p / R_rcs) = 10 N m per rad/s
% whatever the inertia. On the hopper (Ixx 12.7) that's a 0.8 rad/s loop;
% on the CTV (Ixx 0.0123) 810 rad/s, far past the props' 0.1 s spin-up,
% and the CTV spun up in roll to hundreds of rad/s. R_rcs is set so the
% roll-rate loop crosses at half the prop bandwidth (1 / (2 tau)).
Ixx = IN.vehicle.mass.dry_I(1, 1);
wb = 0.5 / IN.vehicle.rcs.tau;                     % rad/s, roll-rate loop bandwidth
Qp = 1e6 / 10000;                                  % LinerizationMaster's roll-rate weight
% Thrust: the hopper's thrust weight on a vehicle 55x lighter gives a
% ~750 N/m height gain, a 17 rad/s vertical loop instead of the hopper's
% 2.3, and flying on the estimate every mm of height noise moved thrust
% by 0.75 N, every step: the chattering thrust, sampled by the IMU, left a
% velocity bias that walked the estimate off. Scaling the weight by the
% mass ratio squared keeps K / m, so the CTV gets the hopper's vertical
% dynamics.
m_ratio = vehicleConfig('hopper').mass.wet / IN.vehicle.mass.wet;
IN.gnc.lqr.R_diag = [3e7 * m_ratio^2, 8e6, 8e6, Qp / (Ixx * wb)^2 * 10000] / 10000;   % pitch, yaw: the hopper's
% Roll angle: with the hopper's quaternion weights the angle gain fell with
% R_rcs to nearly nothing, so a steady roll torque (the TVC's, ~0.04 N m
% flying on the estimate) held the CTV turning at ~0.8 rad/s. An angle
% loop at w_th: K_th = Ixx * w_th^2; the quaternion moves half the roll
% angle, so its weight is 4 * R_rcs * K_th^2. A steady torque M then
% leaves a fixed roll offset M / K_th (~13 deg for 0.04 N m), harmless
% on an axisymmetric vehicle.
w_th = 0.8 * wb;                                   % rad/s
IN.gnc.lqr.roll_angle_weight = 4 * IN.gnc.lqr.R_diag(4) * (Ixx * w_th^2)^2;

% ---- Wind: calm. ASSUMPTION: the CTV flies in calm air; the hopper's
% 10 m/s default puts ~11 N of side force (~40 % of weight) on it
IN.wind.uwind = 0;                                 % m/s
IN.wind.vwind = 0;                                 % m/s

% ---- Mission: 4 m hop with limits scaled to the CTV ----
IN.mission.target_altitude = 4;                    % m  [req]
IN.mission.max_ascent_vel = 2;                     % m/s [req] ASSUMPTION (CTV inputs had 26)
IN.mission.max_descent_vel = 2;                    % m/s [req] ASSUMPTION
% The hopper's landing limits before they moved to the GNC PDR values,
% kept so the CTV's requirements don't loosen with the hopper's
IN.mission.max_impact_vel = 0.5;                   % m/s [req] ASSUMPTION
IN.mission.max_landing_offset = 0.5;               % m  [req] ASSUMPTION
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
    if isfield(s, 'mass_factor'),    IN.mass_factor = s.mass_factor; end   % vehicleSetup scales the dry mass
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

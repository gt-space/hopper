function W = navWorkspace(W, variants)
%NAVWORKSPACE  Add what the NAVIGATION choices read to W.
%   W = navWorkspace(W, variants)   W from hopperWorkspace (after
%                                   sensorWorkspace), variants from
%                                   variantConfig
%
%   W.nav_params  .F (nav_frames), .S (the Sensors block's sim_sensor_params:
%                 the sigmas that made the noise go into the buses), .geom
%                 (lidar beams and mounts, GPS antennas: sim body, from the
%                 base), .dt (filter step)
%   W.nav_init    InitConditions bus: x0 (the true pad state as the nav
%                 state) and P0
%   W.nav_Qc, W.nav_dt, W.nav_g, W.nav_reset   the filter's other inputs
%
%   The tuning stays the nav team's: P0, Qc and the filter rate come from
%   their filter_sim_initializer (subsystems/navigation), run as is. To try
%   other tuning, set W.nav_Qc (13 x 13) or W.nav_P0 (19 x 19) before
%   buildSimInput; they are kept, like W.faults. The torque block of Qc
%   comes from the vehicle config (VC.nav.torque_psd), since it scales with
%   the vehicle's inertia. The
%   lidar and GPS geometry comes from the vehicle config instead of their
%   hard-coded frames, so the filter and the sensor model agree.

if ~isfield(W, 'sensor_params')
    error('navWorkspace:sensors', 'Run sensorWorkspace first: the filter reads the shared Sensors block.');
end
N = navInitializer();
S = W.sensor_params.S;
if abs(N.filter_st - W.sensor_params.dt) > 1e-12
    error('navWorkspace:rate', 'The filter step (%g s) and the Sensors block step (%g s, the IMU period) must match.', ...
        N.filter_st, W.sensor_params.dt);
end

P = navParams(W, S, N.filter_st);
F = P.F; dcm_pad = F.R';
W.nav_params = P;

cg0 = [W.cg_init; 0; 0];
x_pad = [sim_pos_offset(cg0) + [0; 0; -W.cg_init]; zeros(6, 1); nav_quat_from_dcm(dcm_pad)];
if ~isfield(W, 'nav_P0'), W.nav_P0 = N.P0; end
if ~isfield(W, 'nav_Qc')
    W.nav_Qc = N.qc;
    VC = W.IN.vehicle;
    if isfield(VC, 'nav') && isfield(VC.nav, 'torque_psd')
        % the vehicle's torque disturbance (it scales with inertia), body
        % axes [roll pitch yaw] into the nav frame
        R = P.F.R;
        W.nav_Qc(4:6, 4:6) = R * diag(VC.nav.torque_psd) * R';
    end
end
W.nav_init = struct('x0', nav_state_from_sim(x_pad, dcm_pad, cg0, F), 'P0', W.nav_P0);
W.nav_dt = N.filter_st;
W.nav_g = W.IN.const.g0;
W.nav_reset = 0;                    % ValueType reset is a double
end

function N = navInitializer()
% the nav team's filter_sim_initializer, run in its own folder (it adds
% its subfolders to the path relative to there)
here = fileparts(which('filter_sim_initializer'));
if isempty(here)
    error('navWorkspace:initializer', 'filter_sim_initializer (subsystems/navigation) is not on the path.');
end
run(fullfile(here, 'filter_sim_initializer.m'));
N.filter_st = filter_st;
N.P0 = P0;
N.qc = qc;
end

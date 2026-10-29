function [y, st] = fdir_frame(in, P, st)
% FDIR_FRAME  One 100 Hz frame: monitors -> executive -> mode manager ->
% touchdown -> abort reference (GNC PDR slide 78 order).
%
%   in  .t  .meas .valid .stamp (sensor_layout: the shared Sensors block,
%       sensor_frame, sampled at this frame)  .x (13-state the flight
%       computer has: truth, or the navigation estimate)  .f_b (3)
%       .dcm_be (3x3)
%       .thrust  .gimbal_fb (2)  .rcs_torque  .mass  .u_nom (4)
%       .u_cmd (4) the command FDIR sent (after command select), last frame
%   P   fdirWorkspace's fdir_params: .S .mon .modes .cat .h_cut .abort
%       .faults .inject .ref0_err .thrust_max .thrust_min .g .cutoff_alt
%   st  state (fdir_frame_state)
%
%   y.abort .safed .degrade .nav_abort  response flags for command select
%   y.x_ref     abort reference (13)
%   y.mode .phase .level .first_fault .no_arm .mask   for logging
%   y.status    fault states (catalog order)
%   y.raw       this frame's monitor results (catalog order), for tuning
%   y.nav_nis   normalized NIS per sensor (baro gyro gps lidar dgps), -1
%               without NAVIGATION or a new sample
%   y.meas      the sensor measurements this frame saw
%
%#codegen

x = in.x(:);
t = in.t;

% ---- sensor ingest (sensor_frame did the sensors and the sensor faults) ----
meas = in.meas(:); valid = logical(in.valid(:)); stamp = in.stamp(:);

% ---- G&C and actuator signals ----
alt = -x(3);
body_x_ned = in.dcm_be' * [1; 0; 0];
k.tilt = acosd(min(1, max(-1, -body_x_ned(3))));     % body x (thrust axis) from vertical, deg
k.alt = alt;
k.rates = x(7:9);
k.vel_ned = x(4:6);
k.thrust_cmd = in.u_cmd(1);          % what FDIR sent last frame, not the nominal command
k.thrust = in.thrust;
k.thrust_max = P.thrust_max;
k.thrust_min = P.thrust_min;
k.gimbal_cmd = in.u_cmd(2:3);
k.gimbal_fb = in.gimbal_fb(:);
k.rcs_cmd = in.u_cmd(4);
k.rcs_torque = in.rcs_torque;
k.ref0_err = P.ref0_err;

mode_now = st.modes.mode;
M = fdir_modes();
phase = M.MODE_PHASE(mode_now);
[raw, metric, st.mon] = fdir_gnc_monitors(meas, valid, stamp, t, phase, k, P.mon, st.mon);
% estimator checks (rows 33-37) on the navigation filter, when it runs
nav_nis = -ones(5, 1);
if P.nav_active
    [rn, mn, nav_nis, st.nav] = fdir_nav_monitors(meas, valid, stamp, in.xn(:), in.P, in.cg(:), P.mon, st.nav);
    raw = raw | rn;
    metric = metric + mn;
end

% ---- executive ----
eu.t = t; eu.phase = phase; eu.alt = alt; eu.tw = in.thrust / max(in.mass * P.g, 1e-6); eu.reset = false;
[ex, st.exec] = fdir_executive(raw, metric, eu, P.cat, P.h_cut, st.exec);

% ---- touchdown and mode manager ----
thrust_cut = alt <= P.cutoff_alt && x(6) > 0;
L = sensor_layout();
lid = sort(meas(L.LIDAR));
td.lidar_median = (lid(2) + lid(3)) / 2;
f_now = norm(in.f_b);
td.accel_impulse = abs(f_now - st.prev_f);
st.prev_f = f_now;
td.vz = -x(6);
td.thrust_cut = thrust_cut;
[touchdown, ~, st.modes] = fdir_touchdown(td, P.modes, mode_now, st.modes);

mu.t = t; mu.alt = alt; mu.vz = -x(6);
mu.thrust_ratio = eu.tw;
mu.response = uint8(max(0, ex.level));
mu.uncontrollable = ex.uncontrollable;
mu.touchdown = touchdown;
mu.thrust_cut = thrust_cut;
mu.no_arm = ex.no_arm;
mu.cmd = struct('init', false, 'arm', false, 'start', false, 'launch', false, ...
    'abort', false, 'terminate', false, 'reset', false);
[mo, st.modes] = fdir_mode_manager(mu, P.modes, st.modes);

% ---- abort reference ----
[x_ref, st.abort] = fdir_abort_reference(x, t, mo.abort, P.abort, st.abort);

y.abort = mo.abort;
y.safed = mo.safed;
y.nav_abort = mo.nav_abort;
y.degrade = ex.degrade && ~mo.abort && ~mo.safed;
y.x_ref = x_ref;
y.mode = double(mo.mode);
y.phase = double(mo.phase);
y.level = double(ex.level);
y.first_fault = ex.first_fault;
y.no_arm = ex.no_arm;
y.mask = ex.mask;
y.status = double(ex.status);
y.raw = double(raw);
y.nav_nis = nav_nis;
y.meas = meas;
end


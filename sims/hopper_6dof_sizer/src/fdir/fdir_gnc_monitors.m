function [raw, metric, m] = fdir_gnc_monitors(meas, valid, stamp, t, phase, k, c, m)
% FDIR_GNC_MONITORS  One frame of the GNC monitors, in fault-table order.
%
%   meas, valid, stamp  sensors after the fault-injection gates
%   t, phase            sim time; mode phase (fdir_modes)
%   k                   G&C and actuator signals this frame:
%       .alt (m)  .tilt (deg) body axis from vertical  .rates (3, rad/s, truth)
%       .vel_ned (3)  .thrust_cmd .thrust .thrust_max .thrust_min (N)
%       .gimbal_cmd (2) .gimbal_fb (2) (rad)  .rcs_cmd .rcs_torque (N m)
%       .ref0_err (m) t = 0 reference vs pad state
%   c                   thresholds (gnc_monitor_params)
%   m                   state (gnc_monitor_state), returned updated
%
%   raw     37 x 1 logical: this frame's check failed (catalog row order)
%   metric  37 x 1 statistic for CUSUM rows (0 elsewhere)
%
% Rows 33-37 (estimator checks) are fdir_nav_monitors', set by fdir_frame
% when NAVIGATION is on; false here.
% Checks from FDIR CoDR slides 20-28 and GNC PDR slides 80-83.
%
%#codegen

L = sensor_layout();
N = 37;
raw = false(N, 1);
metric = zeros(N, 1);
on_pad = phase <= uint8(2);
g = 9.80665;

acc = meas(L.ACCEL); gyr = meas(L.GYRO);
imu_nan = any(isnan(acc)) || any(isnan(gyr));

% ---- 1 IMU loss ----
raw(1) = imu_nan || ~valid(L.S_IMU);

% ---- 2-3 IMU saturation ----
if ~imu_nan
    raw(2) = any(abs(acc) > c.imu_fs_frac * c.accel_fs);
    raw(3) = any(abs(gyr) > c.imu_fs_frac * c.gyro_fs);
end

% ---- 4-5 IMU stuck: windowed variance below the noise floor ----
if stamp(L.S_IMU) ~= m.prev_stamp(L.S_IMU) && ~imu_nan
    m.imu_hist = [m.imu_hist(:, 2:end), [acc; gyr]];
    m.imu_n = min(m.imu_n + 1, size(m.imu_hist, 2));
end
if m.imu_n >= c.stuck_window
    v = var(m.imu_hist, 0, 2);
    raw(4) = any(v(1:3) < c.accel_var_floor);
    raw(5) = any(v(4:6) < c.gyro_var_floor);
end

% ---- 6 IMU rate of change ----
prev_acc = m.prev_meas(L.ACCEL);
if ~imu_nan && ~any(isnan(prev_acc))
    raw(6) = any(abs(acc - prev_acc) > c.accel_jump);
end

% ---- 7 IMU timing: dt outside [0.5, 1.5] x nominal ----
age = t - stamp(L.S_IMU);
raw(7) = age > 1.5 * c.imu_period + 1e-9;
if stamp(L.S_IMU) ~= m.prev_stamp(L.S_IMU) && m.prev_stamp(L.S_IMU) >= 0
    dt = stamp(L.S_IMU) - m.prev_stamp(L.S_IMU);
    raw(7) = raw(7) || dt < 0.5 * c.imu_period - 1e-9 || dt > 1.5 * c.imu_period + 1e-9;
end

% ---- 8-12 Lidar per beam (blind zone handled by the phase mask) ----
% Only inside the rated envelope: lidar is used near the ground (PDR 80),
% so no return at altitude is expected, not a fault.
in_range = k.alt / cosd(c.lidar_cant) < c.lidar_envelope * c.lidar_r_max;
all_bad = in_range;
for b = 1:4
    r = meas(L.LIDAR(b));
    bad = isnan(r) || ~valid(L.S_LIDAR(b)) || r < c.lidar_r_min || r > c.lidar_r_max;
    rp = m.prev_meas(L.LIDAR(b));
    dts = stamp(L.S_LIDAR(b)) - m.prev_stamp(L.S_LIDAR(b));
    if ~bad && ~isnan(rp) && dts > 0
        bad = abs(r - rp) / dts > c.lidar_rate;
    end
    raw(7 + b) = bad && in_range;
    all_bad = all_bad && bad;
end
raw(12) = all_bad;

% ---- 13-15 GPS ----
raw(13) = meas(L.GPS_FIX) < 3 || meas(L.GPS_SATS) < c.gps_min_sats || ~valid(L.S_GPS);
raw(14) = meas(L.GPS_AGE) > c.gps_max_age;
pg = meas(L.GPS_POS);
if stamp(L.S_GPS) ~= m.prev_stamp(L.S_GPS) && ~any(isnan(pg))
    if m.gps_prev_t >= 0 && ~any(isnan(m.gps_prev_pos))
        dtg = stamp(L.S_GPS) - m.gps_prev_t;
        jump = abs((pg - m.gps_prev_pos) - k.vel_ned(:) * dtg);
        raw(15) = any(jump > c.gps_jump_k * sqrt(2) * c.gps_sigma);
    end
    m.gps_prev_pos = pg;
    m.gps_prev_t = stamp(L.S_GPS);
end
raw(15) = raw(15) || any(isnan(pg));

% ---- 16 Barometer: range, rate ----
hb = meas(L.BARO);
raw(16) = isnan(hb) || ~valid(L.S_BARO) || hb < c.baro_lo || hb > c.baro_hi;
hbp = m.prev_meas(L.BARO);
dtb = stamp(L.S_BARO) - m.prev_stamp(L.S_BARO);
if ~raw(16) && ~isnan(hbp) && dtb > 0
    raw(16) = abs(hb - hbp) / dtb > c.baro_rate;
end

% ---- 17 Magnetometer ----
bm = meas(L.MAG);
raw(17) = any(isnan(bm)) || ~valid(L.S_MAG) || norm(bm) < c.mag_lo || norm(bm) > c.mag_hi;

% ---- 18-22 Guidance and control ----
raw(18) = k.tilt > c.tilt_e1;
raw(19) = k.tilt > c.tilt_e2 || any(abs(k.rates(2:3)) > c.rate_w2);   % roll rate is fault 22 (R2)
% tilt low-passed before differencing: an estimated attitude (NAVIGATION
% gmekf) carries noise that a raw 100 Hz difference turns into deg/s
if isnan(m.prev_tilt_err), m.prev_tilt_err = k.tilt; m.tilt_f = k.tilt; m.tilt_slow = k.tilt; end
tf = m.tilt_f + (c.frame / c.limit_cycle_tau) * (k.tilt - m.tilt_f);
er = (tf - m.tilt_f) / c.frame;
m.tilt_f = tf;
% oscillation amplitude: tilt about its short-term average, over the
% window. A limit cycle swings degrees; flying on an estimate jitters
% tenths (22 sign changes in 2 s but < 0.5 deg about a 0.2 s average)
m.tilt_slow = m.tilt_slow + (c.frame / c.limit_cycle_avg) * (k.tilt - m.tilt_slow);
m.tilt_osc = [m.tilt_osc(2:end), abs(tf - m.tilt_slow)];
sgn = 0;
if er > c.limit_cycle_deadband, sgn = 1; elseif er < -c.limit_cycle_deadband, sgn = -1; end
m.err_rate_sign = [m.err_rate_sign(2:end), sgn];
changes = 0; last = 0;
W = numel(m.err_rate_sign);
for j = W - c.limit_cycle_window + 1:W
    s = m.err_rate_sign(j);
    if s ~= 0
        if last ~= 0 && s ~= last, changes = changes + 1; end
        last = s;
    end
end
raw(20) = changes > c.limit_cycle_count && ...
    max(m.tilt_osc(end - c.limit_cycle_window + 1:end)) > c.limit_cycle_amp;
m.prev_tilt_err = k.tilt;
raw(21) = k.thrust_cmd > c.throttle_frac * k.thrust_max;
raw(22) = abs(k.rates(1)) > c.roll_rate_max;
raw(23) = on_pad && k.ref0_err > c.ref0_tol;

% ---- 24-29 Actuators ----
at_lim = abs(k.gimbal_fb(:)) >= c.gimbal_at_limit * c.gimbal_max;
m.gimbal_hist = [m.gimbal_hist(:, 2:end), at_lim];
Wg = size(m.gimbal_hist, 2);
duty1 = sum(m.gimbal_hist(:, Wg - c.gimbal_duty_win1 + 1:Wg), 2) / c.gimbal_duty_win1;
duty2 = sum(m.gimbal_hist(:, Wg - c.gimbal_duty_win2 + 1:Wg), 2) / c.gimbal_duty_win2;
raw(24) = any(duty1 > c.gimbal_duty1);
raw(25) = any(duty2 > c.gimbal_duty2);
res = abs(k.gimbal_cmd(:) - k.gimbal_fb(:));
raw(26) = any(res > c.gimbal_residual);
raw(27) = any(at_lim & (res > c.gimbal_hardover * 2 * c.gimbal_max));
% expected roll torque through EDF_Dynamics, in its order: first-order lag
% on the command, quantize, torque rate limit, clamp at max_thrust (the
% block clamps the torque at the thrust value, so the clamp rarely acts)
tq_max = c.rcs_max_thrust * c.rcs_arm;
m.rcs_filtered = m.rcs_filtered + (c.frame / c.rcs_tau) * (k.rcs_cmd - m.rcs_filtered);
tq_q = c.rcs_resolution * round(m.rcs_filtered / c.rcs_resolution);
m.rcs_limited = m.rcs_limited + min(max(tq_q - m.rcs_limited, -c.rcs_rate * c.frame), c.rcs_rate * c.frame);
tq_exp = min(max(m.rcs_limited, -c.rcs_max_thrust), c.rcs_max_thrust);
raw(28) = abs(tq_exp - k.rcs_torque) > c.rcs_residual_frac * tq_max;
% expected thrust: the command within the throttle range (the engine can't
% go below its minimum), through the throttle
T_cmd = min(max(k.thrust_cmd, k.thrust_min), k.thrust_max);
% valve rate limit (on the pad and at startup it follows the engine)
if on_pad
    m.thrust_expected = k.thrust;
else
    m.thrust_expected = m.thrust_expected + min(max(T_cmd - m.thrust_expected, -c.thrust_rate * c.frame), c.thrust_rate * c.frame);
end
raw(29) = abs(k.thrust - m.thrust_expected) > c.thrust_residual_frac * k.thrust_max;

% ---- 30-32 IBIT on the pad ----
if on_pad && ~imu_nan
    raw(30) = abs(norm(acc) - g) > c.pad_g_tol || any(abs(gyr) > c.pad_gyro_max);
end
if on_pad
    if any(isnan(m.pad_lidar0)), m.pad_lidar0 = meas(L.LIDAR); end
    raw(31) = any(abs(meas(L.LIDAR) - m.pad_lidar0) > c.pad_lidar_tol);
    raw(32) = any(abs(pg) > c.pad_gps_tol);
end

m.prev_meas = meas;
m.prev_stamp = stamp;
end

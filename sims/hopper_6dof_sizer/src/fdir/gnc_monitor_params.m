function c = gnc_monitor_params(IN, S)
% GNC_MONITOR_PARAMS  Thresholds of the GNC monitors (fdir_gnc_monitors).
%
%   c = gnc_monitor_params(IN, S)   S = sim_sensor_params(IN)
%
% Provisional values: the GNC PDR (slide 82, "Open") sets the thresholds
% and persistence times from Monte Carlo false-alarm vs missed-detection
% statistics. IN.fdir.monitors overrides any field.

c.frame = S.frame;

% ---- IMU (CoDR 20) ----
c.imu_fs_frac    = 0.95;                       % saturation at 0.95 full scale
c.accel_fs       = S.imu.accel_fs;
c.gyro_fs        = S.imu.gyro_fs;
c.stuck_window   = 20;                         % samples for the windowed variance
c.accel_var_floor = (0.1 * min(S.imu.accel_sigma))^2;
c.gyro_var_floor  = (0.1 * min(S.imu.gyro_sigma))^2;
c.accel_jump     = 20;                         % m/s^2 per sample (physical jerk bound)
c.imu_period     = max(S.imu.period, S.frame);   % s: FDIR sees one IMU sample per frame at most (IMU 500 Hz, frame 100 Hz)

% ---- Lidar per beam (CoDR 21) ----
c.lidar_r_min = S.lidar.r_min;
c.lidar_r_max = S.lidar.r_max;
c.lidar_rate  = 10 / cosd(S.lidar.cant + 15) + 2;   % m/s: v_z,max / cos(theta_max) + margin
c.lidar_cant  = S.lidar.cant;
c.lidar_envelope = 0.9;                       % checked below 90 % of the rated range

% ---- GPS (CoDR 22) ----
c.gps_min_sats = 7;                            % > 6
c.gps_max_age  = 0.2;                          % s
c.gps_jump_k   = 5;
c.gps_sigma    = S.gps.sigma(:);
c.gps_period   = S.gps.period;

% ---- Barometer / mag (PDR 80) ----
c.baro_lo = -500; c.baro_hi = 9000;
c.baro_rate = 50;                              % m/s
c.mag_lo = 20; c.mag_hi = 80;                  % uT field magnitude

% ---- Guidance and control (PDR 82, CoDR 27-28) ----
c.tilt_e1       = 10;                          % deg, loss-of-control bound
c.tilt_e2       = 30;                          % deg, attitude envelope
c.rate_w2       = deg2rad(60);                 % rad/s, attitude envelope (pitch and yaw rates)
c.roll_rate_max = deg2rad(90);                 % rad/s
c.throttle_frac = 0.9;                         % thrust margin
c.limit_cycle_window = 200;                    % frames (2 s)
c.limit_cycle_count  = 12;                     % sign changes in the window
c.limit_cycle_deadband = 0.5;                  % deg/s
c.limit_cycle_amp = 0.5;                       % deg, oscillation about the short-term average a limit cycle must reach
c.limit_cycle_avg = 0.2;                       % s, that average
c.limit_cycle_tau = 0.05;                      % s, tilt low-pass before the rate (passes the few-Hz cycles it looks for)

% ---- Actuators (CoDR 28, PDR 83) ----
c.gimbal_max       = IN.tvc.max_gimbal_angle;
c.gimbal_at_limit  = 0.98;                     % fraction of max counted as "at the limit"
c.gimbal_duty_win1 = 100; c.gimbal_duty1 = 0.3;    % 1 s, high duty -> R1
c.gimbal_duty_win2 = 200; c.gimbal_duty2 = 0.8;    % 2 s, sustained -> R2
c.gimbal_residual  = deg2rad(3);               % stuck / slow
c.gimbal_hardover  = 0.5;                      % residual beyond half the range
c.rcs_residual_frac = 0.3;                     % of max roll torque
c.rcs_max_thrust   = IN.rcs.max_thrust;       % N (EDF_Dynamics model)
c.rcs_arm          = IN.vehicle.rcs.arm;      % m, EDF_Dynamics moment arm
c.rcs_rate         = 25;                      % N m/s, EDF_Dynamics rate limiter (on the torque)
c.rcs_tau          = IN.rcs_time_constant;    % s
c.rcs_resolution   = IN.rcs.resolution;       % N m, EDF_Dynamics quantizer
c.thrust_rate    = IN.throttle_valve.rate_limit;   % N/s, throttle valve rate limit
c.thrust_residual_frac = 0.1;                  % stuck thrust: |T - T_cmd| > 10 % of max

% ---- IBIT on the pad (CoDR 20-22) ----
c.pad_g_tol    = 0.5;                          % m/s^2, |f| ~ g
c.pad_gyro_max = 0.05;                         % rad/s
c.pad_lidar_tol = 0.05;                        % m
c.pad_gps_tol  = 0.5;                          % m from the NED origin
c.ref0_tol     = 0.05;                         % m, t = 0 reference vs pad state

% ---- Estimator checks, rows 33-37 (fdir_nav_monitors; NAVIGATION on) ----
% 99.99 % chi-square gates on the filter's reported NIS (the filter gates
% outliers at 99 %; FDIR looks for persistent faults, and 99 % on five
% sensor groups at 100 Hz trips the N-of-M by chance)
c.gate_baro  = 15.137; c.gate_gyro = 21.108; c.gate_gps = 27.856;
c.gate_lidar = 15.137; c.gate_dgps = 21.108;
c.lidar_beams_fail = 2;                        % beams failing together before the lidar group counts (one beam is row 9-12)
c.pos_var_max   = 10^2;                       % m^2: position no longer known to 10 m
c.vel_var_max   = 5^2;                        % (m/s)^2
c.att_var_max   = deg2rad(20)^2;              % rad^2 per axis
c.quat_norm_tol = 1e-3;
c.speed_max     = 100;                        % m/s, far beyond the hopper's envelope
c.rate_max      = 5;                          % rad/s
c.shadow_sigmas = 5;                           % row 36: main vs shadow, combined sigmas
c.shadow_dz     = 1.0;                         % m, and at least this far apart in height
c.shadow_dvz    = 0.5;                         % m/s, or in vertical velocity
c.nav_r_scale   = 4;                           % R inflation while the NIS CUSUM (row 35) is latched
c.nav = struct();                             % nav_params, set by fdirWorkspace

if isfield(IN, 'fdir') && isfield(IN.fdir, 'monitors')
    for f = fieldnames(IN.fdir.monitors)'
        c.(f{1}) = IN.fdir.monitors.(f{1});
    end
end
end

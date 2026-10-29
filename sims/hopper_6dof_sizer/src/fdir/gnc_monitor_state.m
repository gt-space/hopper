function m = gnc_monitor_state()
% GNC_MONITOR_STATE  Initial state for fdir_gnc_monitors.
%
%#codegen

L = sensor_layout();
W = 200;                                     % longest window (frames)
m.imu_hist   = zeros(6, 20);                 % last stuck_window IMU samples
m.imu_n      = 0;
m.prev_meas  = NaN(L.N_MEAS, 1);
m.prev_stamp = -ones(L.N_SENSORS, 1);
m.imu_dt_bad = false;
m.gps_prev_pos = NaN(3, 1);
m.gps_prev_t   = -1;
m.gimbal_hist  = false(2, W);                % at-limit flags, pitch / yaw
m.err_rate_sign = zeros(1, W);               % sign of tilt-error rate
m.prev_tilt_err = 0;
m.rcs_filtered  = 0;                         % expected roll torque after the fan lag (N m)
m.rcs_limited   = 0;                         % ... and after the torque rate limit (N m)
m.thrust_expected = 0;                       % thrust command after the throttle rate limit (N)
m.pad_lidar0    = NaN(4, 1);                 % beam ranges at the first pad frame
end

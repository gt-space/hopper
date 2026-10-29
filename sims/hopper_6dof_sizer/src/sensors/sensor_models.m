function [meas, truth, valid, stamp, st] = sensor_models(x, noise, t, S, st)
% SENSOR_MODELS  One frame of the vehicle's sensors (sensor_layout order).
%
%   x      truth: .pos_ned (3, the vehicle's base reference point, z = 0
%          on the pad), .dcm_be (3x3, NED to body), .w_b (3),
%          .f_b (3, specific force, body)
%   noise  20 x 1 standard-normal samples for this frame (seeded Random
%          Number block, so runs repeat exactly)
%   t      sim time (s)
%   S      sim_sensor_params
%   st     state (sensor_state): held samples and their times
%
%   meas   noisy measurements; truth the same channels without noise;
%   valid  device valid flags; stamp time of each sensor's latest sample.
% Each sensor samples at its own period and holds between samples, as the
% flight computer's ring buffers would.
%
%#codegen

L = sensor_layout();
pos = x.pos_ned(:); C = x.dcm_be;

% ---- noise-free channels ----
tr = zeros(L.N_MEAS, 1);
tr(L.ACCEL) = x.f_b(:);
tr(L.GYRO)  = x.w_b(:);
tr(L.BARO)  = -pos(3);
for k = 1:4
    az = S.lidar.azimuth(k);
    d_b = [-cosd(S.lidar.cant); sind(S.lidar.cant) * cosd(az); sind(S.lidar.cant) * sind(az)];   % beam, body (down the thrust axis, canted out)
    r_b = [S.lidar.mount_x; S.lidar.radius * cosd(az); S.lidar.radius * sind(az)];               % mount from the base reference
    p_n = pos + C' * r_b;
    d_n = C' * d_b;
    if d_n(3) > 1e-3
        tr(L.LIDAR(k)) = max(0, -p_n(3)) / d_n(3);   % to the pad plane z = 0 (NED, down +)
    else
        tr(L.LIDAR(k)) = S.lidar.r_max + 1;          % pointing up or level: no return
    end
end
tr(L.GPS_POS)  = pos;
tr(L.GPS_FIX)  = S.gps.fix;
tr(L.GPS_SATS) = S.gps.sats;
tr(L.GPS_AGE)  = 0;
tr(L.MAG)      = C * S.mag.field_ned(:);

% ---- noisy channels ----
sig = zeros(L.N_MEAS, 1);
sig(L.ACCEL) = S.imu.accel_sigma(:);
sig(L.GYRO)  = S.imu.gyro_sigma(:);
sig(L.BARO)  = S.baro.sigma;
sig(L.LIDAR) = S.lidar.sigma;
sig(L.GPS_POS) = S.gps.sigma(:);
sig(L.MAG)   = S.mag.sigma;
noisy = tr + sig .* noise(:);
noisy(L.ACCEL) = min(max(noisy(L.ACCEL), -S.imu.accel_fs), S.imu.accel_fs);
noisy(L.GYRO)  = min(max(noisy(L.GYRO),  -S.imu.gyro_fs),  S.imu.gyro_fs);
noisy(L.LIDAR) = max(noisy(L.LIDAR), 0);

% ---- sample and hold, per sensor ----
period = [S.imu.period; S.baro.period; repmat(S.lidar.period, 4, 1); S.gps.period; S.mag.period];
for s = 1:L.N_SENSORS
    if t - st.stamp(s) >= period(s) - 1e-9 || st.stamp(s) < 0
        ch = L.CHANNEL_SENSOR(:) == s;
        st.meas(ch)  = noisy(ch);
        st.truth(ch) = tr(ch);
        st.stamp(s)   = t;
    end
end
meas = st.meas;
truth = st.truth;
meas(L.GPS_AGE) = t - st.stamp(L.S_GPS);
truth(L.GPS_AGE) = meas(L.GPS_AGE);
stamp = st.stamp;

valid = true(L.N_SENSORS, 1);
for k = 1:4   % lidar returns only inside its rated range
    r = meas(L.LIDAR(k));
    valid(L.S_LIDAR(k)) = r >= S.lidar.r_min && r <= S.lidar.r_max;
end
end

function [meas, truth, valid, stamp, st] = sensor_models(x, noise, t, S, st)
% SENSOR_MODELS  One frame of the vehicle's sensors (sensor_layout order).
%
%   x      truth: .pos_ned (3, the controller's x(1:3): the CG position
%          shifted to the base as if upright, sim_pos_offset; z = 0 on the
%          pad), .vel_ned (3, CG velocity), .dcm_be (3x3, NED to body),
%          .w_b (3), .f_b (3, specific force, body), .cg (3, Mass_Properties
%          CG from the base, body)
%   noise  26 x 1 standard-normal samples for this frame (seeded Random
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
cg = x.cg(:); w = x.w_b(:);
p_cg = pos - sim_pos_offset(cg);            % CG, NED
% a point fixed in the body (r_b from the base, body axes): position, velocity
at = @(r_b) p_cg + C' * (r_b - cg);
vat = @(r_b) x.vel_ned(:) + C' * cross(w, r_b - cg);

% ---- noise-free channels ----
tr = zeros(L.N_MEAS, 1);
tr(L.ACCEL) = x.f_b(:);
tr(L.GYRO)  = x.w_b(:);
tr(L.BARO)  = -pos(3);
[d_lid, r_lid] = lidar_geometry(S.lidar);
for k = 1:4
    p_n = at(r_lid(:, k));
    d_n = C' * d_lid(:, k);
    if d_n(3) > 1e-3
        tr(L.LIDAR(k)) = max(0, -p_n(3)) / d_n(3);   % to the pad plane z = 0 (NED, down +)
    else
        tr(L.LIDAR(k)) = S.lidar.r_max + 1;          % pointing up or level: no return
    end
end
ant = S.gps.antenna;                         % 3 x 2, from the base, body
tr(L.GPS_POS)  = at(ant(:, 1));
tr(L.GPS_VEL)  = vat(ant(:, 1));
tr(L.DGPS)     = C' * (ant(:, 2) - ant(:, 1));
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
sig(L.GPS_VEL) = S.gps.vel_sigma(:);
sig(L.DGPS)  = S.dgps.sigma(:);
noisy = tr + sig .* noise(:);
noisy(L.ACCEL) = min(max(noisy(L.ACCEL), -S.imu.accel_fs), S.imu.accel_fs);
noisy(L.GYRO)  = min(max(noisy(L.GYRO),  -S.imu.gyro_fs),  S.imu.gyro_fs);
noisy(L.LIDAR) = max(noisy(L.LIDAR), 0);

% ---- sample and hold, per sensor ----
period = [S.imu.period; S.baro.period; repmat(S.lidar.period, 4, 1); S.gps.period; S.mag.period; S.dgps.period];
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

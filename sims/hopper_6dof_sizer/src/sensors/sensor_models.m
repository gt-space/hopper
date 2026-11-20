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
noise = noise(:);
pos = x.pos_ned(:); C = x.dcm_be;
cg = x.cg(:); w = x.w_b(:);
persistent nav_frame_R;
if isempty(nav_frame_R)
    dcm_pad = angle2dcm(0, pi/2, 0);
    nav_frame_R = nav_frames(dcm_pad, 0).R;
end
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

% Build the state and frames expected by the navigation Src_sim models.
navState = zeros(20, 1);
navState(1:3) = p_cg;
navState(4:6) = x.vel_ned(:);
navState(7:9) = nav_frame_R * w;
C_nav = nav_frame_R * C;
navState(10:13) = nav_quat_from_dcm(C_nav);
navParams = nav_sensor_params();
accelerationTrue = C' * x.f_b(:) + [0; 0; navParams.imu.gravity_scalar];

baroState = navState;
baroState(3) = pos(3);

lidarFrames.heightCOM = 0;
lidarFrames.A1 = completeFrame(nav_frame_R * d_lid(:, 1));
lidarFrames.A2 = completeFrame(nav_frame_R * d_lid(:, 2));
lidarFrames.A3 = completeFrame(nav_frame_R * d_lid(:, 3));
lidarFrames.A4 = completeFrame(nav_frame_R * d_lid(:, 4));
lidarFrames.r1 = nav_frame_R * (r_lid(:, 1) - cg);
lidarFrames.r2 = nav_frame_R * (r_lid(:, 2) - cg);
lidarFrames.r3 = nav_frame_R * (r_lid(:, 3) - cg);
lidarFrames.r4 = nav_frame_R * (r_lid(:, 4) - cg);
gpsFrame.r1 = nav_frame_R * (ant(:, 1) - cg);
gpsFrame.r2 = nav_frame_R * (ant(:, 2) - cg);

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
noisy(L.GYRO)  = min(max(noisy(L.GYRO),  -S.imu.gyro_fs), S.imu.gyro_fs);
noisy(L.LIDAR) = max(noisy(L.LIDAR), 0);

navInjectedNoise = sig .* noise(:);
navInjectedNoise(L.ACCEL) = nav_frame_R * navInjectedNoise(L.ACCEL);
navInjectedNoise(L.GYRO) = nav_frame_R * navInjectedNoise(L.GYRO);
navInjectedNoise(L.MAG) = nav_frame_R * navInjectedNoise(L.MAG);

% ---- sample and hold, per sensor ----
period = [S.imu.period; S.baro.period; repmat(S.lidar.period, 4, 1); S.gps.period; S.mag.period; S.dgps.period];
for s = 1:L.N_SENSORS
    if t - st.stamp(s) >= period(s) - 1e-9 || st.stamp(s) < 0
        switch s
            case L.S_IMU
                [accelBus, angvelBus] = imu_model(navState, accelerationTrue, ...
                    zeros(3, 1), zeros(3, 1), t, navInjectedNoise);
                noisy(L.ACCEL) = nav_frame_R' * accelBus.accel;
                noisy(L.ACCEL) = min(max(noisy(L.ACCEL), -S.imu.accel_fs), S.imu.accel_fs);
                noisy(L.GYRO) = nav_frame_R' * angvelBus.angvel;
                noisy(L.GYRO) = min(max(noisy(L.GYRO), -S.imu.gyro_fs), S.imu.gyro_fs);

            case L.S_BARO
                baroBus = barometer_model(baroState, 0, t, sig.*noise);
                noisy(L.BARO) = baroBus.height;

            case L.S_LIDAR(1)
                lidarBus = lidar_system_model(navState, t, lidarFrames, sig.*noise);
                noisy(L.LIDAR) = max(lidarBus.distance, 0);
                for k = 1:4
                    d_n = C' * d_lid(:, k);
                    if d_n(3) <= 1e-3
                        noisy(L.LIDAR(k)) = max(tr(L.LIDAR(k)) + sig(L.LIDAR(k))*noise(L.LIDAR(k)), 0);
                    end
                end

            case L.S_GPS
                gpsBus = gps_model(navState, t, gpsFrame, sig.*noise);
                noisy(L.GPS_POS) = gpsBus.pos;
                noisy(L.GPS_VEL) = gpsBus.vel;

            case L.S_MAG
                magBus = mag_model(navState, t, navInjectedNoise, 10 * S.mag.field_ned(:));
                noisy(L.MAG) = 0.1 * (nav_frame_R' * magBus.bfield_meas);

            case L.S_DGPS
                dgpsBus = diffGps_model(navState, t, gpsFrame, sig.*noise);
                noisy(L.DGPS) = dgpsBus.basisVec;
        end
        ch = L.CHANNEL_SENSOR(:) == s;
        st.meas(ch)  = noisy(ch);
        st.truth(ch) = tr(ch);
        st.stamp(s)  = t;
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

function A = completeFrame(direction)
d = direction(:) / norm(direction);
a = [1; 0; 0];
if abs(d(1)) > 0.9, a = [0; 1; 0]; end
e2 = cross(d, a);
e2 = e2 / norm(e2);
A = [d, e2, cross(d, e2)];
end

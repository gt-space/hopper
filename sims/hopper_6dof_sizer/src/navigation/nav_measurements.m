function [accel, angvel, baro, gps, lidar, dgps, mag, lidarFrames, gpsFrame, massProps, torque] = ...
    nav_measurements(meas, valid, stamp, cg, mass, I_b, M_b, P, mask, r_scale)
%NAV_MEASUREMENTS  One sample of the shared sensor model as the filter's buses.
%   [accel, angvel, baro, gps, lidar, dgps, mag, lidarFrames, gpsFrame, ...
%    massProps, torque] = nav_measurements(meas, valid, stamp, cg, mass, I_b, M_b, P, mask, r_scale)
%     meas, valid, stamp  sensor_frame output (sensor_layout)
%     cg      Mass_Properties CG [x; y; z] from the base, sim body (m)
%     mass    current mass (kg); I_b inertia about the CG, sim body (kg m^2)
%     M_b     control torque, sim body (N m)
%     P       nav_params (navWorkspace): .F (nav_frames), .S
%             (sim_sensor_params: the sigmas that made the noise), .geom
%             (lidar beams and mounts, GPS antennas: sim body, from the base)
%     mask    FDIR's sensor mask (fdir_targets bits, fdir_frame y.nav_mask):
%             a masked sensor goes to the filter as invalid, so it is not
%             fused (imu, imu_accel, imu_gyro, lidar1-4, lidar, gps, baro,
%             mag; gps covers the differential GPS too)
%     r_scale FDIR's measurement-noise inflation (1: none), on every
%             covariance
%   Outputs are structs with the fields of NavigationDictionary.sldd's
%   buses (accelBus, angvelBus, baroBus, gpsBus, lidarBus, diffGpsBus,
%   magBus, lidarFrameBus, gpsFrameBus, massProperties) and the torque.
%
%   Everything goes into the nav frame (nav_frames): body vectors are
%   rotated by F.R, positions move to the CG and the nav origin (the CG's
%   rest point on the pad, F.heightCOM above it). Each bus carries the
%   covariance of the noise sensor_models added, so the filter's
%   measurement noise matches the sim's.
%
%   The barometer reads the base's height above the pad, the filter
%   expects the CG's height above its rest point: they differ by the CG's
%   travel as propellant burns (slow, absorbed by the filter's baro bias
%   state) and by cg*(1 - cos tilt).
%
%#codegen

L = sensor_layout();
F = P.F; S = P.S; G = P.geom;
on = @(k) bitand(uint32(mask), bitshift(uint32(1), k - 1)) == 0;   % target k not masked
T = targetBits();
R = F.R;
r_cg = cg(:);
up = [0; 0; F.heightCOM];   % nav origin offset: p_nav = p_ned + up

accel.accel      = R * meas(L.ACCEL);
accel.covariance = r_scale * R * diag(S.imu.accel_sigma(:).^2) * R';
accel.timestamp  = stamp(L.S_IMU);
accel.validity   = logical(valid(L.S_IMU)) && on(T.imu) && on(T.imu_accel);

angvel.angvel     = R * meas(L.GYRO);
angvel.covariance = r_scale * R * diag(S.imu.gyro_sigma(:).^2) * R';
angvel.timestamp  = stamp(L.S_IMU);
angvel.validity   = logical(valid(L.S_IMU)) && on(T.imu) && on(T.imu_gyro);

baro.height     = meas(L.BARO);
baro.covariance = r_scale * S.baro.sigma^2;
baro.timestamp  = stamp(L.S_BARO);
baro.validity   = logical(valid(L.S_BARO)) && on(T.baro);

gps.pos         = meas(L.GPS_POS) + up;   % antenna 1, NED
gps.pcovariance = r_scale * diag(S.gps.sigma(:).^2);
gps.vel         = meas(L.GPS_VEL);
gps.vcovariance = r_scale * diag(S.gps.vel_sigma(:).^2);
gps.timestamp   = stamp(L.S_GPS);
gps.validity    = logical(valid(L.S_GPS)) && on(T.gps);

lidar.distance   = meas(L.LIDAR);
lidar.covariance = r_scale * repmat(S.lidar.sigma^2, 4, 1);
lidar.timestamp  = max(stamp(L.S_LIDAR));
lidar.validity   = logical(valid(L.S_LIDAR(:))) & meas(L.LIDAR) <= S.lidar.r_max & ...
    [on(T.lidar1); on(T.lidar2); on(T.lidar3); on(T.lidar4)] & on(T.lidar);

dgps.basisVec   = meas(L.DGPS);           % antenna 1 -> 2, NED
dgps.covariance = r_scale * diag(S.dgps.sigma(:).^2);
dgps.timestamp  = stamp(L.S_DGPS);
dgps.validity   = logical(valid(L.S_DGPS)) && on(T.gps);

mag.bfield_meas = R * meas(L.MAG);
mag.covariance  = r_scale * R * diag(repmat(S.mag.sigma^2, 3, 1)) * R';
mag.bfield_E    = S.mag.field_ned(:);         % Earth field, NED (mag_update_model)
mag.timestamp   = stamp(L.S_MAG);
mag.validity    = logical(valid(L.S_MAG)) && on(T.mag);

% Lidar k: beam direction in A(:,1) (lidar_update_model), mount from the CG
A = zeros(3, 3, 4); r = zeros(3, 4);
for k = 1:4
    d = R * G.lidar_dir(:, k);
    A(:, :, k) = completeFrame(d);
    r(:, k) = R * (G.lidar_mount(:, k) - r_cg);
end
lidarFrames.A1 = A(:, :, 1); lidarFrames.A2 = A(:, :, 2);
lidarFrames.A3 = A(:, :, 3); lidarFrames.A4 = A(:, :, 4);
lidarFrames.r1 = r(:, 1); lidarFrames.r2 = r(:, 2);
lidarFrames.r3 = r(:, 3); lidarFrames.r4 = r(:, 4);
lidarFrames.heightCOM = F.heightCOM;

gpsFrame.r1 = R * (G.gps_antenna(:, 1) - r_cg);
gpsFrame.r2 = R * (G.gps_antenna(:, 2) - r_cg);

massProps.m = mass;
massProps.J = R * I_b * R';
torque = R * M_b(:);
end

function A = completeFrame(d)
% orthonormal frame with d as its first column (only that column is used)
d = d / norm(d);
a = [1; 0; 0];
if abs(d(1)) > 0.9, a = [0; 1; 0]; end
e2 = cross(d, a); e2 = e2 / norm(e2);
A = [d, e2, cross(d, e2)];
end

function T = targetBits()
% fdir_targets positions of the sensor targets
T = struct('imu', 1, 'imu_accel', 2, 'imu_gyro', 3, 'lidar1', 4, 'lidar2', 5, 'lidar3', 6, 'lidar4', 7, ...
           'lidar', 8, 'gps', 9, 'baro', 10, 'mag', 11);
end

function [accel, angvel, baro, gps, lidar, dgps, mag, lidarFrames, gpsFrame, massProps, torque] = ...
    nav_measurements(meas, valid, stamp, cg, mass, I_b, M_b, P)
%NAV_MEASUREMENTS  One sample of the shared sensor model as the filter's buses.
%   [accel, angvel, baro, gps, lidar, dgps, mag, lidarFrames, gpsFrame, ...
%    massProps, torque] = nav_measurements(meas, valid, stamp, cg, mass, I_b, M_b, P)
%     meas, valid, stamp  sensor_frame output (sensor_layout)
%     cg      Mass_Properties CG [x; y; z] from the base, sim body (m)
%     mass    current mass (kg); I_b inertia about the CG, sim body (kg m^2)
%     M_b     control torque, sim body (N m)
%     P       nav_params (navWorkspace): .F (nav_frames), .S
%             (sim_sensor_params: the sigmas that made the noise), .geom
%             (lidar beams and mounts, GPS antennas: sim body, from the base)
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
R = F.R;
r_cg = cg(:);
up = [0; 0; F.heightCOM];   % nav origin offset: p_nav = p_ned + up

accel.accel      = R * meas(L.ACCEL);
accel.covariance = R * diag(S.imu.accel_sigma(:).^2) * R';
accel.timestamp  = stamp(L.S_IMU);
accel.validity   = logical(valid(L.S_IMU));

angvel.angvel     = R * meas(L.GYRO);
angvel.covariance = R * diag(S.imu.gyro_sigma(:).^2) * R';
angvel.timestamp  = stamp(L.S_IMU);
angvel.validity   = logical(valid(L.S_IMU));

baro.height     = meas(L.BARO);
baro.covariance = S.baro.sigma^2;
baro.timestamp  = stamp(L.S_BARO);
baro.validity   = logical(valid(L.S_BARO));

gps.pos         = meas(L.GPS_POS) + up;   % antenna 1, NED
gps.pcovariance = diag(S.gps.sigma(:).^2);
gps.vel         = meas(L.GPS_VEL);
gps.vcovariance = diag(S.gps.vel_sigma(:).^2);
gps.timestamp   = stamp(L.S_GPS);
gps.validity    = logical(valid(L.S_GPS));

lidar.distance   = meas(L.LIDAR);
lidar.covariance = repmat(S.lidar.sigma^2, 4, 1);
lidar.timestamp  = max(stamp(L.S_LIDAR));
lidar.validity   = logical(valid(L.S_LIDAR(:))) & meas(L.LIDAR) <= S.lidar.r_max;

dgps.basisVec   = meas(L.DGPS);           % antenna 1 -> 2, NED
dgps.covariance = diag(S.dgps.sigma(:).^2);
dgps.timestamp  = stamp(L.S_DGPS);
dgps.validity   = logical(valid(L.S_DGPS));

mag.bfield_meas = R * meas(L.MAG);
mag.covariance  = R * diag(repmat(S.mag.sigma^2, 3, 1)) * R';
mag.bfield_E    = norm(S.mag.field_ned);
mag.timestamp   = stamp(L.S_MAG);
mag.validity    = logical(valid(L.S_MAG));

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

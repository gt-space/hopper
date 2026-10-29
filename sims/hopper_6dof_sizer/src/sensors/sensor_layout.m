function L = sensor_layout()
% SENSOR_LAYOUT  Index map of the Sensors block's output vectors.
%
% The block puts out three vectors per frame:
%   meas   26 x 1  measurements (noisy, or noise-free for the truth copy)
%   valid   9 x 1  device valid flag per sensor
%   stamp   9 x 1  time of the latest sample per sensor (s)
% FDIR, the fault-injection gates and navigation (nav_measurements) index them
% through this map, so the layout lives in one place.
%
%#codegen

% ---- meas ----
L.ACCEL    = 1:3;     % IMU specific force, body axes (m/s^2)
L.GYRO     = 4:6;     % IMU body rates (rad/s)
L.BARO     = 7;       % barometric altitude above the pad (m)
L.LIDAR    = 8:11;    % lidar beam ranges, beams 1-4 (m)
L.GPS_POS  = 12:14;   % GPS position of antenna 1, NED from the pad origin (m)
L.GPS_FIX  = 15;      % fix type (3 = 3D)
L.GPS_SATS = 16;      % satellites used
L.GPS_AGE  = 17;      % solution age (s)
L.MAG      = 18:20;   % magnetic field, body axes (uT)
L.GPS_VEL  = 21:23;   % GPS velocity of antenna 1, NED (m/s)
L.DGPS     = 24:26;   % differential GPS baseline, antenna 1 -> 2, NED (m)
L.N_MEAS   = 26;

% ---- valid / stamp: one entry per sensor ----
L.S_IMU   = 1;
L.S_BARO  = 2;
L.S_LIDAR = 3:6;
L.S_GPS   = 7;
L.S_MAG   = 8;
L.S_DGPS  = 9;
L.N_SENSORS = 9;

% sensor of each meas channel
L.CHANNEL_SENSOR = uint8([1 1 1 1 1 1 2 3 4 5 6 7 7 7 7 7 7 8 8 8 7 7 7 9 9 9]);
end

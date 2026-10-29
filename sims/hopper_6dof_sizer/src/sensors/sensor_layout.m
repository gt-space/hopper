function L = sensor_layout()
% SENSOR_LAYOUT  Index map of the Sensors block's output vectors.
%
% The block puts out three vectors per frame:
%   meas   20 x 1  measurements (noisy, or noise-free for the truth copy)
%   valid   8 x 1  device valid flag per sensor
%   stamp   8 x 1  time of the latest sample per sensor (s)
% FDIR, the fault-injection gates and (later) navigation all index them
% through this map, so the layout lives in one place.
%
%#codegen

% ---- meas ----
L.ACCEL    = 1:3;     % IMU specific force, body axes (m/s^2)
L.GYRO     = 4:6;     % IMU body rates (rad/s)
L.BARO     = 7;       % barometric altitude above the pad (m)
L.LIDAR    = 8:11;    % lidar beam ranges, beams 1-4 (m)
L.GPS_POS  = 12:14;   % GPS position, NED from the pad origin (m)
L.GPS_FIX  = 15;      % fix type (3 = 3D)
L.GPS_SATS = 16;      % satellites used
L.GPS_AGE  = 17;      % solution age (s)
L.MAG      = 18:20;   % magnetic field, body axes (uT)
L.N_MEAS   = 20;

% ---- valid / stamp: one entry per sensor ----
L.S_IMU   = 1;
L.S_BARO  = 2;
L.S_LIDAR = 3:6;
L.S_GPS   = 7;
L.S_MAG   = 8;
L.N_SENSORS = 8;

% sensor of each meas channel
L.CHANNEL_SENSOR = uint8([1 1 1 1 1 1 2 3 4 5 6 7 7 7 7 7 7 8 8 8]);
end

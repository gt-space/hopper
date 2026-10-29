function F = fault_codes()
% FAULT_CODES  Targets and modes of the fault table (faultScenario).
%
% The fault table is an N x 6 numeric matrix so it can be a Simulink
% parameter and a Monte Carlo variable:
%   [target, mode, t_on, t_off, magnitude, param]
% faultScenario builds it from these names.
%
%#codegen

% ---- targets ----
F.T_ACCEL    = 1;   % IMU accelerometer, all axes (param: axis 1-3, 0 = all)
F.T_GYRO     = 2;   % IMU gyro (param: axis)
F.T_BARO     = 3;
F.T_LIDAR    = 4;   % param: beam 1-4, 0 = all
F.T_GPS      = 5;
F.T_MAG      = 6;
F.T_GIMBAL_P = 11;  % TVC pitch
F.T_GIMBAL_Y = 12;  % TVC yaw
F.T_THROTTLE = 13;
F.T_RCS      = 14;

% ---- modes ----
F.M_LOSS       = 1;    % NaN, device invalid
F.M_BIAS       = 2;    % + magnitude
F.M_DRIFT      = 3;    % + magnitude * (t - t_on)
F.M_SATURATION = 4;    % value pinned at magnitude
F.M_STUCK      = 5;    % hold the value at onset
F.M_NOISE      = 6;    % + magnitude-amplitude burst
F.M_STALE      = 7;    % timestamp and value stop updating
F.M_POOR_FIX   = 8;    % GPS: fix type magnitude, satellites param
F.M_SLOW       = 9;    % actuator command rate-limited to magnitude per s
F.M_HARDOVER   = 10;   % actuator command pinned at magnitude
F.M_DEAD       = 11;   % actuator output zero
F.M_REDUCED    = 12;   % actuator output scaled by magnitude

F.NCOL = 6;
end

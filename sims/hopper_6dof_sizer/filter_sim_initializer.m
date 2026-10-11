%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter_sim_initializer.m
% 
% Sets all parameters for the navigation filter simulation, including the
% sensor parameters.
% 
% Author: David Meisinger
% Date: 10/14/26
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Filter parameters
filter_st = 0.002;      % Filter sample time (seconds)

%% Sensor parameters

% IMU parameters
imu_st = 1/960;         % sample time (s)

% GPS parameters
gps_st = 0.1;           % sample time (s)

% Barometer parameters
barometer_st = 0.05;    % sample time (s)

% LIDAR parameters
lidar_st = 0.01;        % sample time (s)

% LIDAR 1, 2, 3 mounting parameters

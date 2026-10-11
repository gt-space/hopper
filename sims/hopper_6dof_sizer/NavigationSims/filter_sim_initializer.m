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

%% Add relevant file paths
addpath("Filter Sims\");
addpath("Src\");
addpath("Filters\");

%% Filter parameters
filter_st = 1/480;      % Filter sample time (seconds)
% NOTE: This comes from IMU sample rate.
% IT IS ASSUMED THE FILTER ITERATES (PROPAGATES) AT THE IMU SENSE RATE.

% Initial conditions
% state: xx, xy, xz, vx, vy, vz, P, Q, R, q0, q1, q2, q3, ax, ay, az, gx,
% gy, gz, b
x0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% TODO! TUNE THIS
P0 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
% Note covariance is 19x19, not 20x20, because covariance uses a 3-element
% rotation representation

%% Sensor sampling and noise parameters

% IMU parameters
imu_st = filter_st;     % sample time (s)

% GPS parameters
gps_st = 0.1;           % sample time (s)

% Barometer parameters
barometer_st = 0.05;    % sample time (s)

% LIDAR parameters
lidar_st = 0.01;        % sample time (s)

%% Sensor frame parameters

% IMU Frame Transformation


% LIDAR 1, 2, 3 mounting parameters

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

% TODO! TUNE THIS COVARIANCE
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
% --- IMU Sensor Parameters ---
imu_params.accel_bias = [0.01; -0.01; 0.02];     % Example accelerometer bias (m/s^2)
imu_params.gyro_bias = [0.001; 0.002; -0.001];   % Example gyroscope bias (rad/s)
imu_params.accel_noise_std = [0.05; 0.05; 0.05]; % Accelerometer noise standard deviation
imu_params.gyro_noise_std = [0.005; 0.005; 0.005];% Gyroscope noise standard deviation
imu_params.gravity_body = [0; 0; -9.81];         % Gravity component in body frame
imu_params.packet_invalid_prob = 0.0;            % Set to 0 for now as requested

% GPS parameters
gps_st = 0.1;           % sample time (s)
% --- GPS Sensor Parameters (u-blox ZED-F9P) ---
gps_params.r_b = [0.1; 0.0; -0.2];             % Antenna lever arm from CG (3x1) (m)
gps_params.bias = [0.005; -0.005; 0.01];       % Slowly drifting error (b_k) (m)[cite: 1]
gps_params.multipath = [0.002; 0.001; -0.003]; % Pad bounce error (b_mp,k) (m)[cite: 1]
gps_params.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter standard deviation (eta_p) (m)[cite: 1]
gps_params.packet_invalid_prob = 0.0;          % Set to 0 for now

% Barometer parameters
barometer_st = 0.05;    % sample time (s)
% --- Barometer Sensor Parameters ---
baro_params.bias = 0.0;                 % Height bias (m)
baro_params.noise_std = 0.1;            % Noise standard deviation (m)
baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now

% LIDAR parameters
lidar_st = 0.01;        % sample time (s)

%% Sensor frame parameters

% IMU Frame Transformation
% LIDAR 1, 2, 3 mounting parameters
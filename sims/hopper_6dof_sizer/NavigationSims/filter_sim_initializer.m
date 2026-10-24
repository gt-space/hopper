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
addpath("Src_filter\");
addpath("Filters\");

%% Filter parameters
filter_st = 1/480;      % Filter sample time (seconds)
% NOTE: This comes from IMU sample rate.
% IT IS ASSUMED THE FILTER ITERATES (PROPAGATES) AT THE IMU SENSE RATE.
% Initial conditions

% state: xx, xy, xz, vx, vy, vz, P, Q, R, q0, q1, q2, q3, ax, ay, az, gx,
% gy, gz, b
x0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

%% GMEKF continuous-time process-noise spectral densities

% Translational acceleration disturbance: 3 channels
accelerationNoisePSD = imu_params.accel_noise_std.^2 * imu_st;

% requires calibration
torqueNoisePSD = zeros(3,1);
accelBiasNoisePSD = zeros(3,1);
gyroBiasNoisePSD = zeros(3,1);
baroBiasNoisePSD = 0;

Qc = blkdiag( ...
    diag(accelerationNoisePSD), ...
    diag(torqueNoisePSD), ...
    diag(accelBiasNoisePSD), ...
    diag(gyroBiasNoisePSD), ...
    baroBiasNoisePSD);

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
imu_params.accel_noise_std = [0.05; 0.05; 0.05]; % Accelerometer noise standard deviation
imu_params.gyro_noise_std = [0.005; 0.005; 0.005];% Gyroscope noise standard deviation
imu_params.gravity_scalar = 9.81;                % Gravity acceleration
imu_params.packet_invalid_prob = 0.0;            % Set to 0 for now as requested

% GPS parameters
gps_st = 0.1;           % sample time (s)
% --- GPS Sensor Parameters (u-blox ZED-F9P) ---
gps_params.r_b = [0.1; 0.0; -0.2];             % Antenna lever arm from CG (3x1) (m)
gps_params.multipath = [0.002; 0.001; -0.003]; % Pad bounce error (b_mp,k) (m)[cite: 1]
gps_params.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter standard deviation (eta_p) (m)[cite: 1]
gps_params.packet_invalid_prob = 0.0;          % Set to 0 for now

% Barometer parameters
barometer_st = 0.05;    % sample time (s)
% --- Barometer Sensor Parameters ---
baro_params.noise_std = 0.1;            % Noise standard deviation (m)
baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now

% LIDAR parameters
lidar_st = 0.01;        % sample time (s)

%% Sensor frame parameters

% IMU Frame Transformation


% LIDAR 1, 2, 3, 4 mounting frame information
% TODO: get these params from CAD accurately
lidar_frames.A1 = lidarRotMat(0, 20);
lidar_frames.A2 = lidarRotMat(90, 20);
lidar_frames.A3 = lidarRotMat(180, 20);
lidar_frames.A4 = lidarRotMat(270, 20);
lidar_frames.r1 = [0.25 0 0]';
lidar_frames.r2 = [0 0.25 0]';
lidar_frames.r3 = [-0.25 0 0]';
lidar_frames.r4 = [0 -0.25 0]';

%% Function to build rotation matrix for each LIDAR unit vector
% phi: angle from x about z-axis of body frame (deg)
% theta: angle between lidar sensor direction and body z-axis (deg)
function A_Li = lidarRotMat(phi, theta)
psi = theta - 90;
cp = cosd(phi);
sp = sind(phi);
ct = cosd(psi);
st = sind(psi);
A3 = [cp -sp 0;
      sp  cp 0;
      0   0  1];
A2 = [ ct 0 st;
       0  1 0;
      -st 0 ct];
A_Li = A3*A2;
end
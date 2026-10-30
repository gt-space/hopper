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
addpath("Src_sim\");
addpath("Src_filter\");
addpath("Filters\");

%% Filter parameters
filter_st = 1/500;      % Filter sample time (seconds)
% NOTE: This comes from IMU sample rate. ADIS16500 output rate is
% 2000/(DEC_RATE+1) Hz, so 500 Hz needs DEC_RATE = 3 in the FSW driver.
% IT IS ASSUMED THE FILTER ITERATES (PROPAGATES) AT THE IMU SENSE RATE.
% Initial conditions

% state: xx, xy, xz, vx, vy, vz, P, Q, R, q0, q1, q2, q3, ax, ay, az, gx,
% gy, gz, b
x0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

% TODO! TUNE THIS COVARIANCE
P0 = diag([ ...
    0.1*ones(1,3), ...      % pos (m^2)
    0.01*ones(1,3), ...     % vel (m/s)^2
    zeros(1,3), ...         % omega: no longer a state
    (deg2rad(3))^2*ones(1,3), ...   % attitude, ~3 deg
    (0.05)^2*ones(1,3), ... % accel bias, ~0.05 m/s^2
    (deg2rad(1))^2*ones(1,3), ...   % gyro bias, ~1 deg/s
    0.01 ]);                % baro bias

% Note covariance is 19x19, not 20x20, because covariance uses a 3-element
% rotation representation

% Enforce some stuff for flying the gyro
P0(7:9,:) = 0;
P0(:,7:9) = 0;

%% Sensor sampling and noise parameters

% IMU parameters
imu_st = filter_st;     % sample time (s)
% --- IMU Sensor Parameters (ADIS16500, data sheet Rev. A, Table 1) ---
% Per-sample white noise at the output rate: std = density * sqrt(rate/2)
imu_params_init.accel_noise_density = [880e-6; 880e-6; 732e-6];   % (m/s^2)/sqrt(Hz), x y z
imu_params_init.gyro_noise_density = deg2rad([6.1e-3; 6.1e-3; 7.0e-3]); % (rad/s)/sqrt(Hz), x y z
imu_params_init.accel_bias_stability = [125e-6; 125e-6; 134e-6];  % In-run bias stability (m/s^2)
imu_params_init.gyro_bias_stability = deg2rad([7.5; 8.1; 4.9])/3600; % In-run bias stability (rad/s)
imu_params_init.accel_noise_std = imu_params_init.accel_noise_density * sqrt(1/(2*imu_st)); % Accelerometer noise standard deviation (m/s^2)
imu_params_init.acovariance = diag(imu_params_init.accel_noise_std.^2);
imu_params_init.gyro_noise_std = imu_params_init.gyro_noise_density * sqrt(1/(2*imu_st)); % Gyroscope noise standard deviation (rad/s)
imu_params_init.gcovariance = diag(imu_params_init.gyro_noise_std.^2);
imu_params_init.gravity_scalar = 9.81;                % Gravity acceleration
imu_params_init.packet_invalid_prob = 0.0;            % Set to 0 for now as requested

% GPS parameters
gps_st = 0.1;           % sample time (s)
% --- GPS Sensor Parameters (u-blox ZED-F9P) ---
% PLACEHOLDER antenna spacing: 0.5 m separation per David; tune with CAD
% later. Keep the baseline horizontal (body x-y, the N-E plane when
% vertical) so DGPS observes roll about the vertical axis.
gps_params_init.r1 = [0.25; 0.0; -0.2];             % Antenna 1 vector from COM (3x1) (m)
gps_params_init.r2 = [-0.25; 0.0; -0.2];            % Antenna 2 vector from COM (3x1) (m)
% gps_params.bias = [0; 0; 0];                   % Constant position bias (m)
% gps_params.multipath = [0.002; 0.001; -0.003]; % Pad bounce error (b_mp,k) (m)[cite: 1]
% gps_params.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter standard deviation (eta_p) (m)[cite: 1]
% gps_params.vel_noise_std = [0.05, 0.05, 0.01]; % Velocity jitter standard deviation
% gps_params.pos_covariance = diag([0.02^2; 0.02^2; 0.05^2]); % Filter measurement covariance
% gps_params.vel_covariance = diag([0.05^2; 0.05^2; 0.1^2]); % Filter measurement covariance
% gps_params.packet_invalid_prob = 0.0;          % Set to 0 for now
% % Provisional GPS velocity noise standard deviations (m/s).
% % Replace with receiver specifications.
% gps_params.vel_noise_std = [0.1; 0.1; 0.1];

% Differential GPS parameters
diffGps_st = 0.1;       % sample time (s), PLACEHOLDER: same as gps_st until receiver is chosen
% --- Differential GPS Sensor Parameters (antenna 1 -> antenna 2 baseline in NED) ---
% Antenna placement comes from gps_params.r1 / gps_params.r2 (gps_frame).
% PLACEHOLDER noise: carrier-phase baseline error is roughly constant in
% metres, so heading error ~ noise_std/baseline. 7 mm horizontal on the
% 0.5 m baseline is ~0.8 deg. Replace with receiver specs and CAD spacing.
% diffGps_params.noise_std = [0.007; 0.007; 0.012];          % Baseline jitter per NED axis (m)
% diffGps_params.covariance = diag(diffGps_params.noise_std.^2); % Filter measurement covariance
% diffGps_params.packet_invalid_prob = 0.0;                  % Set to 0 for now

% Barometer parameters
barometer_st = 0.05;    % sample time (s)
% --- Barometer Sensor Parameters ---
% baro_params.noise_std = 0.1;            % Noise standard deviation (m)
% baro_params.covariance = 0.1^2;         % Filter measurement covariance
% baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now
% baro_params.bias = 0.5;                 % Placeholder

% LIDAR parameters
lidar_st = 0.01;        % sample time (s)
% --- Barometer Sensor Parameters (SJ-GU-TF-Luna A03)---
% lidar_params.noise_std_1 = 0.06; % Distance measurement jitter (cm) for 0.2-3m
% lidar_params.noise_std_2 = 0.02; % Distance measurement jitter (proportion distance) for 3-8m
% lidar_params.packet_invalid_prob = 0.0; % set to zero for now
% Note- any distances over 8m will be not valid it appears

%% GMEKF continuous-time process-noise spectral densities

% Accel noise and the bias errors are expressed in the inertial frame in the
% filter, so use the largest body-axis value on all three channels.

% Translational acceleration disturbance: 3 channels
accelerationNoisePSD = max(imu_params_init.accel_noise_std)^2 * imu_st * ones(3,1);
imu_params_init.gyro_R = diag(imu_params_init.gyro_noise_std.^2);

% Bias random walks from the in-run bias stability. TUNING STARTING POINT:
% the data sheet gives no rate random walk, so assume the bias wanders by
% about its in-run stability over biasCorrelationTime.
biasCorrelationTime = 100;  % (s)
accelBiasNoisePSD = max(imu_params_init.accel_bias_stability)^2 / biasCorrelationTime * ones(3,1);
gyroBiasNoisePSD = max(imu_params_init.gyro_bias_stability)^2 / biasCorrelationTime * ones(3,1);
gyroNoisePSD = imu_params_init.gyro_noise_std.^2 * imu_st;

% requires calibration
baroBiasNoisePSD = 0.001;

qc = blkdiag( ...
    diag(accelerationNoisePSD), ...
    diag(gyroNoisePSD), ...
    diag(accelBiasNoisePSD), ...
    diag(gyroBiasNoisePSD), ...
    baroBiasNoisePSD);

%% Sensor frame parameters

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
lidar_frames.heightCOM = 1;

% GPS mounting parameters
gps_frame.r1 = gps_params_init.r1;
gps_frame.r2 = gps_params_init.r2;

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
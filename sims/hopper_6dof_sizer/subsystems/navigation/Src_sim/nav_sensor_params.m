%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nav_sensor_params.m
%
% One place for the sensor models' parameters (Src_sim). Each model reads
% its struct once into its persistent copy; a test can override them.
%
%   P = nav_sensor_params()          current parameters (defaults below)
%   nav_sensor_params(P)             replace them (e.g. noise off), then
%                                    clear the models' persistents:
%                                    clear gps_model imu_model ...
%   nav_sensor_params('reset')       back to the defaults
%
% P.gps P.imu P.baro P.lidar P.diffGps: the structs the models used to
% build inside themselves, same values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = nav_sensor_params(set)
    persistent current;
    if nargin > 0
        if ischar(set)
            current = [];       % 'reset'
        else
            current = set;
        end
    end
    if isempty(current)
        current = defaults();
    end
    P = current;
end

function P = defaults()
    % --- GPS (u-blox ZED-F9P) ---
    P.gps.r1 = [0.25; 0.0; -0.2];             % Antenna 1 vector from COM (3x1) (m)
    P.gps.r2 = [-0.25; 0.0; -0.2];            % Antenna 2 vector from COM (3x1) (m)
    P.gps.bias = [0; 0; 0];                   % Constant position bias (m)
    P.gps.multipath = [0.002; 0.001; -0.003]; % Pad bounce error (b_mp,k) (m)
    P.gps.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter standard deviation (eta_p) (m)
    P.gps.pos_covariance = diag([0.02^2; 0.02^2; 0.05^2]); % Filter measurement covariance
    P.gps.vel_covariance = diag([0.05^2; 0.05^2; 0.1^2]);  % Filter measurement covariance
    P.gps.packet_invalid_prob = 0.0;          % Set to 0 for now
    % Provisional GPS velocity noise standard deviations (m/s).
    % Replace with receiver specifications.
    P.gps.vel_noise_std = [0.1; 0.1; 0.1];

    % --- IMU (ADIS16500, data sheet Rev. A, Table 1) ---
    imu_st = 1/500;
    P.imu.accel_noise_density = [880e-6; 880e-6; 732e-6];   % (m/s^2)/sqrt(Hz), x y z
    P.imu.gyro_noise_density = deg2rad([6.1e-3; 6.1e-3; 7.0e-3]); % (rad/s)/sqrt(Hz), x y z
    P.imu.accel_bias_stability = [125e-6; 125e-6; 134e-6];  % In-run bias stability (m/s^2)
    P.imu.gyro_bias_stability = deg2rad([7.5; 8.1; 4.9])/3600; % In-run bias stability (rad/s)
    P.imu.accel_noise_std = [880e-6; 880e-6; 732e-6] * sqrt(1/(2*imu_st)); % Accelerometer noise standard deviation (m/s^2)
    P.imu.acovariance = diag([880e-6; 880e-6; 732e-6].^2 * (1/(2*imu_st)));
    P.imu.gyro_noise_std = deg2rad([6.1e-3; 6.1e-3; 7.0e-3]) * sqrt(1/(2*imu_st)); % Gyroscope noise standard deviation (rad/s)
    P.imu.gcovariance = diag(deg2rad([6.1e-3; 6.1e-3; 7.0e-3].^2) * (1/(2*imu_st)));
    P.imu.gravity_scalar = 9.81;                % Gravity acceleration
    P.imu.packet_invalid_prob = 0.0;            % Set to 0 for now as requested

    % --- Barometer ---
    P.baro.noise_std = 0.1;            % Noise standard deviation (m)
    P.baro.covariance = 0.1^2;         % Filter measurement covariance
    P.baro.packet_invalid_prob = 0.0;  % Set to 0 for now
    P.baro.bias = 0.5;                 % Placeholder

    % --- LIDAR (SJ-GU-TF-Luna A03) ---
    P.lidar.noise_std_1 = 0.06;        % Distance measurement jitter (m) for 0.2-3m
    P.lidar.noise_std_2 = 0.02;        % Distance measurement jitter (proportion distance) for 3-8m
    P.lidar.packet_invalid_prob = 0.0; % set to zero for now

    % --- Differential GPS (antenna 1 -> antenna 2 baseline in NED) ---
    P.diffGps.noise_std = [0.007; 0.007; 0.012];          % Baseline jitter per NED axis (m)
    P.diffGps.covariance = diag([0.007; 0.007; 0.012].^2); % Filter measurement covariance
    P.diffGps.packet_invalid_prob = 0.0;                  % Set to 0 for now
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imu_model.m
%
% Function to output the IMU measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
%
% Author: Aditiya Babbar
% Date: 10/11/26
%
% Inputs:
%   x: state vector (20x1)
%   a: acceleration vector (3x1)
%   accelBias: accelerometer measurement bias in body frame (m/s^2)
%   gyroBias: gyroscope measurement bias in body frame (rad/s)
%   simtime: time elapsed (s)
%
% Outputs:
%   accelBus: struct defined in the Navigation Dictionary
%   angvelBus: struct defined in the Navigation Dictionary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accelBus, angvelBus] = imu_model(x, a, accelBias, gyroBias, simtime)

% Build empty output structs from the dictionary definition
accelBus = Simulink.Bus.createMATLABStruct('accelBus');
angvelBus = Simulink.Bus.createMATLABStruct('angvelBus');

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent imu_params;
if isempty(imu_params)
    try
        imu_params = evalin('base', 'imu_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        imu_params.accel_noise_std = [0; 0; 0];
        imu_params.gyro_noise_std = [0; 0; 0];
        imu_params.acovariance = [0.05; 0.05; 0.05].^2*eye(3);
        imu_params.gcovariance = [0.005; 0.005; 0.005].^2*eye(3);
        imu_params.gravity_scalar = 9.81; % Example body frame gravity component
        imu_params.packet_invalid_prob = 0.0;
    end
end

% 1. Forward Timestamp
timestamp = simtime;

% 2. Accelerometer Measurement: A(q)*(a_true + ba + na) - g_vec
% Adjust gravity subtraction based on your coordinate convention/body frame definition
g_vec = [0; 0; -g]; % accel reading in NED frame when IMU is still relative to ground

accel_ideal = a - g_radial;
accel_noise = imu_params.accel_noise_std .* randn(3, 1);
accel = accel_ideal + accelBias + accel_noise;

% 3. Gyroscope Measurement: w_true + bg + ng
% NOTE: Update indices (e.g., x(10:12)) to match your specific 20x1 state vector layout for angular velocity
omega_true = x(7:9);
gyro_noise = imu_params.gyro_noise_std .* randn(3, 1);
angvel = omega_true + gyroBias + gyro_noise;

% 4. Validity Flag Assignment
if rand() < imu_params.packet_invalid_prob
    validity = false;
else
    validity = true;
end

% Build the output structs
accelBus.accel = accel;
accelBus.covariance = imu_params.acovariance;
accelBus.timestamp = timestamp;
accelBus.validity = validity;

angvelBus.angvel = anvel;
angvelBus.covariance = imu_params.gcovariance;
angvelBus.timestamp = timestamp;
angvelBus.validity = validity;

end
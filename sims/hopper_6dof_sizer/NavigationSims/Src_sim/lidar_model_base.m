%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_model_base.m
% 
% Function to output the distance measurement from each LIDAR given the
% true distance between the sensor and ground.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
% 
% Author: David Meisinger
% Date: 10/24/26
% 
% Inputs:
%   dist_true: true distances to measure (m) (4x1)
% 
% Outputs:
%   dist_meas: measured distances with noise (m) (4x1)
%   covariance: variance of each beam for the filter to use (4x1)
%   validity: false if sensor misbehaves
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist_meas, covariance, validity] = lidar_model_base(dist_true)

% Retrieve parameters from the workspace (set by filter_sim_initializer)
% persistent lidar_params;
% if isempty(lidar_params)
%     try
%         lidar_params = evalin('base', 'lidar_params');
%     catch
%         % Fallback defaults if initializer hasn't populated workspace yet
%         lidar_params.noise_std_1 = 0.06; % Distance measurement jitter (cm) for 0.2-3m
%         lidar_params.noise_std_2 = 0.02; % Distance measurement jitter (proportion distance) for 3-8m
%         lidar_params.packet_invalid_prob = 0.0; % set to zero for now
%     end
% end

% 4. Validity Flag Assignment
validity = rand(4,1) > lidar_params.packet_invalid_prob;

% Handle noise and validity according to the true distance of each beam
num_sensors = size(dist_true,1);
noise_std = lidar_params.noise_std_1*ones(num_sensors,1);   % 0.2-3m
far = dist_true > 3;
noise_std(far) = lidar_params.noise_std_2*dist_true(far);   % 3-8m
inRange = dist_true >= 0.2 & dist_true <= 8;
validity = validity & inRange;

dist_noise = noise_std.*randn(num_sensors,1);

% Compute the measured distance
dist_meas = dist_true + dist_noise;

% Reported variance uses the measured distance, as the adaptor would
reported_std = lidar_params.noise_std_1*ones(num_sensors,1);
reported_std(far) = lidar_params.noise_std_2*dist_meas(far);
covariance = reported_std.^2;

end
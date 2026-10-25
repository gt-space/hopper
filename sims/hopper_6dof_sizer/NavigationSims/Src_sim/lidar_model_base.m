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
%   covariance: covariance for the filter to use (4x4)
%   validity: false if sensor misbehaves
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist_meas, covariance, validity] = lidar_model_base(dist_true)

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent lidar_params;
if isempty(lidar_params)
    try
        lidar_params = evalin('base', 'lidar_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        lidar_params.noise_std_1 = 0.06; % Distance measurement jitter (cm) for 0.2-3m
        lider_params.noise_std_2 = 0.02; % Distance measurement jitter (% distance) for 3-8m
        lidar_params.packet_invalid_prob = 0.0; % set to zero for now
    end
end

% 4. Validity Flag Assignment
validity = rand(4,1) > gps_params.packet_invalid_prob;

% Handle noise according to the true distance
if dist_true < 0.2
    validity = zeros(4,1);
end

% For each distance add gaussian noise and determine validity
dist_meas(idx,1) = dist_true(idx,1);
validity(idx,1) = true;

end
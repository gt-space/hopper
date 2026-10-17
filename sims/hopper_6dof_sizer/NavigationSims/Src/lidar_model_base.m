%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_model_base.m
% 
% Function to output the distance measurement from each LIDAR given the
% true distance between the sensor and ground.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
% 
% Author:
% Date: 
% 
% Inputs:
%   dist_true: true distance computed using sensor and vehicle frame (m)
% 
% Outputs:
%   dist_meas: measured distance (m)
%   validity: false if sensor misbehaves
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist_meas = lidar_model_base(dist_true)
    % TODO: IMPLEMENT NOISE AND STUFF
    dist_meas = dist_true;
end
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
%   dist_true: true distances to measure (m) (4x1)
% 
% Outputs:
%   dist_meas: measured distances with noise (m) (4x1)
%   validity: false if sensor misbehaves
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist_meas, validity] = lidar_model_base(dist_true)

    % TODO: IMPLEMENT NOISE AND STUFF
    for idx = 1:size(dist_true,1)
        dist_meas(idx,1) = dist_true(idx,1);
        validity(idx,1) = true;
    end

end
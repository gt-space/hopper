%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar1_model.m
% 
% Function to output the distance measurement from LIDAR 1 accounting for
% vehicle state and the sensor frame. Assume the ground is flat for now.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
% 
% Author:
% Date: 
% 
% Inputs:
%   x: state vector (20x1)
%   simtime: time elapsed (s)
% 
% Outputs:
%   dist: distance measured (m)
%   timestamp: time elapsed (s)
%   validity: boolean (true if valid measurement)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist, timestamp, validity] = lidar1_model(x, time)

    % First calculate the true distance using frame info

    % Calculate the measured distance from the true distance
    dist = lidar_model_base(dist_true);

end
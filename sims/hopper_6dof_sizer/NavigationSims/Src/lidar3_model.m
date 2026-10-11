%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar3_model.m
% 
% Function to output the distance measurement from LIDAR 3 accounting for
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
function [dist, timestamp, validity] = lidar3_model(x, time)

end
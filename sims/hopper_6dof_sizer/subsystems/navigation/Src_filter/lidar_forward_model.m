%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_system_model.m
% 
% Function to predict measurement of four LIDAR sensors given the estimated
% vehicle state and sensor frames.
% Depends on filter_sim_initializer setting the sensor frames.
% 
% Author: David Meisinger
% Date: 10/17/26
% 
% Inputs:
%   x: filter state estimate (20x1)
% 
% Outputs:
%   d_vec: LIDAR sensor distance predictions (4x1) (m)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_vec] = lidar_forward_model(x, time, lidar_frames)

% Load the LIDAR frames


end
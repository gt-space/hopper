%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imu_model.m
% 
% Function to output the IMU measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
% 
% Author:
% Date: 
% 
% Inputs:
%   x: state vector (20x1)
%   a: acceleration vector (3x1)
%   simtime: time elapsed (s)
% 
% Outputs:
%   accel: acceleration vector (3x1) (m)
%   angvel: angular velocity vector (3x1) (radians/s)
%   timestamp: time elapsed (s)
%   validity: boolean (true if valid measurement)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accel, angvel, timestamp, validity] = imu_model(x, a, simtime)

end
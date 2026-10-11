%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gps_model.m
% 
% Function to output the GPS measurement data with noise and bias applied.
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
%   pos: position vector (3x1) (m)
%   timestamp: time elapsed (s)
%   validity: boolean (true if valid measurement)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos, timestamp, validity] = gps_model(x, simtime)

end
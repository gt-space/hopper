%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% barometer_model.m
% 
% Function to output the baro measurement data with noise and bias applied.
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
%   height: height from ground (m)
%   timestamp: time elapsed (s)
%   validity: boolean (true if valid measurement)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [height, timestamp, validity] = barometer_model(x, time)

end
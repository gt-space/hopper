%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar4_model.m
% 
% Function to output the magnetic field measurement from the magnetometer.
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
%   bfield_meas: magnetic field vector measured (3x1)
%   bfield_E: Earth magnetic field vector in NED from bfield model (3x1)
%   timestamp: time elapsed (s)
%   validity: boolean (true if valid measurement)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bfield_meas, bfield_E, timestamp, validity] = MAG_model(x, time)

end
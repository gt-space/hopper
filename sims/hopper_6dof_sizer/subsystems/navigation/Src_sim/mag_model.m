%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar4_model.m
% 
% Function to output the magnetic field measurement from the magnetometer.
% Currently only applies random noise, as our intended use is to find our
% hard-iron calibration offset and apply it in the sensor bridge before
% incorporation into the filter. Soft-iron bias is not simulated.
% 
% Author: David Meisinger
% Date: 11/15/26
% 
% Inputs:
%   x: state vector (20x1)
%   simtime: time elapsed (s)
% 
% Outputs:
%   magBus: struct defined in the Navigation Dictionary
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function magBus = mag_model(x, time)

% Retrieve simulation parameters
persistent mag_params;
if isempty(mag_params)
    P = nav_sensor_params();   % defaults there; tests can override
    mag_params = P.mag;
end

% Extract attitude DCM
q_BI = x(10:13);
A_IB = quaternion_to_rotation(q_BI)';

% Compute measured magnetic field vector with noise applied in body frame
bfield_m = A_IB*mag_params.bfield_E + mag_params.noise_std.*randn(3, 1);

% Build the output struct
magBus.bfield_meas = bfield_m;
magBus.covariance = mag_params.covariance * eye(3);
magBus.bfield_E = mag_params.bfield_E;
magBus.timestamp = time;
magBus.validity = rand() >= mag_params.packet_invalid_prob;

end
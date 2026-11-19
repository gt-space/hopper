%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mag_model.m
% 
% Function to output the magnetic field measurement from the magnetometer.
% Soft-iron and hard-iron calibration terms are not simulated.
% 
% Author: David Meisinger
% Date: 11/15/26
% 
% Inputs:
%   x: state vector (20x1)
%   simtime: time elapsed (s)
%   measurementNoise: optional additive Sensors measurement noise (26x1, uT)
%   fieldNED: optional Earth field override (mGauss, NED)
% 
% Outputs:
%   magBus: struct defined in the Navigation Dictionary (field units: mGauss)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function magBus = mag_model(x, time, measurementNoise, fieldNED)

% Retrieve simulation parameters
persistent mag_params L;
if isempty(mag_params)
    P = nav_sensor_params();   % defaults there; tests can override
    mag_params = P.mag;
    L = sensor_layout();
end

% Extract attitude DCM
q_BI = x(10:13);
A_IB = quaternion_to_rotation(q_BI)';

if nargin >= 4 && ~isempty(fieldNED)
    bfield_E = fieldNED(:);
else
    bfield_E = mag_params.bfield_E(:);
end

hasInjectedNoise = nargin >= 3 && ~isempty(measurementNoise);
if hasInjectedNoise
    magNoise = 10 * measurementNoise(L.MAG);
else
    magNoise = mag_params.noise_std.*randn(3, 1);
end

% Compute measured magnetic field in body axes (mGauss)
bfield_m = A_IB * bfield_E + magNoise;

% Build the output struct
magBus.bfield_meas = bfield_m;
magBus.covariance = mag_params.covariance * eye(3);
magBus.bfield_E = bfield_E;
magBus.timestamp = time;
if hasInjectedNoise
    magBus.validity = true;
else
    magBus.validity = rand() >= mag_params.packet_invalid_prob;
end

end

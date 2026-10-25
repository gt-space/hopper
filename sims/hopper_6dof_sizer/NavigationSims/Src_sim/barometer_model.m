%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% barometer_model.m
%
% Function to output the baro measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
%
% Author: Aditiya Babbar
% Date: 10/11/26
%
% Inputs:
%   x: state vector (20x1)
%   baroBias: barometer measurement bias (m)
%   simtime: time elapsed (s)
%
% Outputs:
%   baroBus: struct defined in the Navigation Dictionary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function baroBus = barometer_model(x, baroBias, time)

% Build an empty output struct from the dictionary definition
baroBus = Simulink.Bus.createMATLABStruct('baroBus');

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent baro_params;
if isempty(baro_params)
    try
        baro_params = evalin('base', 'baro_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        baro_params.noise_std = 0.1;            % Height noise std deviation (m)
        baro_params.covariance = 0.1^2;
        baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now
    end
end

% 1. Forward Timestamp
timestamp = time;

% 2. Compute Height Measurement with Bias and Noise

% Adjust index x(3) to match whatever index represents your true height/altitude in the 20x1 state vector
true_height = x(3);
baro_noise = baro_params.noise_std * randn();
height = true_height + baroBias + baro_noise;

% 3. Validity Flag Assignment
validity = rand() > baro_params.packet_invalid_prob;

% Build the output struct
baroBus.height = height;
baroBus.covariance = baro_params.covariance;
baroBus.timestamp = timestamp;
baroBus.validity = validity;

end
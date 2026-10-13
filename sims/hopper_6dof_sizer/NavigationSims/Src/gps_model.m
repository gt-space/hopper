%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gps_model.m
%
% Function to output the GPS measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
%
% Author: Aditiya Babbar
% Date: 10/11/26
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

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent gps_params;
if isempty(gps_params)
    try
        gps_params = evalin('base', 'gps_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        gps_params.r_b = [0; 0; 0];              % Antenna lever arm (3x1)
        gps_params.bias = [0; 0; 0];             % Slowly drifting error (b_k)[cite: 1]
        gps_params.multipath = [0; 0; 0];        % Pad bounce error (b_mp,k)[cite: 1]
        gps_params.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter (eta_p)[cite: 1]
        gps_params.packet_invalid_prob = 0.0;    % Set to 0 for now
    end
end

% 1. Forward Timestamp
timestamp = simtime;

% 2. Extract States from 20x1 State Vector
p_true = x(1:3);           % True CG position in NED (p^n)[cite: 1]
C_bn = quat2dcm(x(10:14))'; % attitude state

% 3. Compute Position Measurement Equation:
% p_tilde = p^n + C_bn^T * r^b + b_k + b_mp,k + eta_p[cite: 1]
lever_arm_effect = C_bn' * gps_params.r_b;
pos_noise = gps_params.pos_noise_std .* randn(3, 1);
pos = p_true + lever_arm_effect + gps_params.bias + gps_params.multipath + pos_noise;

% 4. Validity Flag Assignment
if rand() < gps_params.packet_invalid_prob
    validity = false;
else
    validity = true;
end
end
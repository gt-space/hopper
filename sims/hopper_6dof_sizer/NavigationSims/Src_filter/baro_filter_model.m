%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baro_filter_model.m
% Computes the predicted barometer measurement and its H matrix row for the filter.
% 
% Inputs:
%   m_X: nominal position vector (3x1)
%   m_bh: nominal barometer bias (scalar)
% 
% Outputs:
%   h_baro: predicted measurement (scalar)
%   H_baro: measurement Jacobian row (1x19)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h_baro, H_baro] = baro_filter_model(m_X, m_bh)
e3 = [0; 0; 1];
% Predicted measurement: z_h = e3^T * X + b_h
h_baro = e3' * m_X + m_bh;
% Measurement Jacobian H for Barometer (1x19 state vector order: [theta, v, X, beta, gamma, lambda])
% H_h = [0(1x3), 0(1x3), -e3^T, 0(1x3), 0(1x3), -1] as noted in the handoff doc[cite: 5]

H_baro = zeros(1, 19);
H_baro(3)  = -1;
H_baro(19) =  1;
end


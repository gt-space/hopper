%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baro_filter_model.m
% Author: Yuma Nagayoshi
% Date: 10/23/2026
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

function [predictedMeasurement, measurementJacobian] = ...
    baro_filter_model(positionEstimate, baroBiasEstimate)

    % Predicted altitude in NED coordinates
    predictedMeasurement = -positionEstimate(3) + baroBiasEstimate;

    % Jacobian with respect to the 19-dimensional error state
    measurementJacobian = zeros(1,19);
    measurementJacobian(3) = -1;
    measurementJacobian(19) = 1;


end
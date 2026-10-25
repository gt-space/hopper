% diffGps_update_model.m
% 
% Function to compute the measurement innovation and Jacobian for the
% filter correct step. This function handle is passed into update.m for
% updating the state estimate.
% 
% Inputs:
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - differential GPS measurement struct
%   params              - Sensor parameters (gpsFrameBus params here)
% 
% Outputs:
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state measurement Jacobian H (1x19)
function [innovation, measurementJacobian] = ...
    diffGps_update_model(nominalState, measurement, params)

end
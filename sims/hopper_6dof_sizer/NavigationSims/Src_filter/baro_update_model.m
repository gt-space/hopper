function [innovation, measurementJacobian] = baro_update_model(nominalState, measurement, params)
% (Yuma's implementation | 10/25/26)
% Barometer forward measurement model for the GMEKF
%
% Inputs:
%
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - Barometer measurement struct
%   params              - Sensor parameters (unused here)
%
% Outputs:
%
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state measurement Jacobian H (1x19)

    positionEstimate = nominalState(1:3);
    baroBiasEstimate = nominalState(20);

    predictedMeasurement = -positionEstimate(3) + baroBiasEstimate; % NED frame

    innovation = predictedMeasurement - measurement.height;
    
    measurementJacobian(3) = -1;
    measurementJacobian(19) = 1;

end
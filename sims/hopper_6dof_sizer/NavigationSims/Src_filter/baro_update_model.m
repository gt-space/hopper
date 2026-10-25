function [innovation, measurementJacobian] = baro_update_model(nominalState, measurement, params)
% (Yuma's implementation | 10/24/26)
% Adapter for the GMEKF barometer measurement update
%
% Inputs:
%
%   nominalState - Nominal state estimate (20x1)
%   measurement  - Barometer measurement struct
%   params       - Sensor parameters (unused here)
%
%
% Outputs:
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state Jacobian (1x19)

    positionEstimate = nominalState(1:3);
    baroBiasEstimate = nominalState(20);
    
    [predictedMeasurement, measurementJacobian] = baro_filter_model(positionEstimate, baroBiasEstimate);
    
    innovation = predictedMeasurement - measurement.height;

end
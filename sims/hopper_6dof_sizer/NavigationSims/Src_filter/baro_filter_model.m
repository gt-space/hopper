function [predictedMeasurement, measurementJacobian] = baro_filter_model(positionEstimate, baroBiasEstimate)
% (Yuma's implementation | 10/24/26)
% Computes predicted barometer measurement and respective measurement Jacobian
%
% Inputs:
%
%   positionEstimate    - Nominal position vector (3x1)
%   baroBiasEstimate    - Barometer Bias Estimate (1x1)
%
% Outputs:
%
%   predictedMeasurement    - Predicted Measurement (1x1)
%   measurementJacobian     - A Measurement Jacobian Row (1x19)


    % Predicted altitude in NED coordinates
    predictedMeasurement = -positionEstimate(3) + baroBiasEstimate;

    % Jacobian with respect to the 19-dimensional error state
    measurementJacobian = zeros(1,19);
    measurementJacobian(3) = -1;
    measurementJacobian(19) = 1;


end
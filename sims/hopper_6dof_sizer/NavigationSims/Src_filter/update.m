function [nominalStatePosterior, errorCovariancePosterior, accepted] = ...
    update(nominalStatePrior, errorCovariancePrior, sensorModel, measurement, params)
% (Yuma's implementation | 10/24/26)
% Generic GMEKF measurement correction
%
% Inputs:
%
%   nominalStatePrior     - Prior nominal state (20x1)
%   errorCovariancePrior  - Prior error covariance (19x19)
%   sensorModel           - Sensor model function handle
%   measurement           - Sensor measurement structure
%   params                - Sensor/update parameters
%
% Outputs:
%
%   nominalStatePosterior    - Corrected nominal state (20x1)
%   errorCovariancePosterior - Corrected covariance (19x19)
%   accepted                 - Boolean flag if measurement passed gating

    [innovation, measurementJacobian] = sensorModel(nominalStatePrior, measurement, params);

    % If sensor models discover a validity issue, they report an empty
    % innovation. This is how we address those.
    if isempty(innovation)
        nominalStatePosterior = nominalStatePrior;
        errorCovariancePosterior = errorCovariancePrior;
        accepted = false;
        return;
    end
    
    measurementNoiseCovariance = measurement.R;
    
    % Innovation gate threshold
    gateThreshold = params.gateThreshold;
    
    % Perform Kalman update
    [errorCorrection, errorCovarianceConditioned, accepted] = ...
        kalman_measurement_update( ...
        errorCovariancePrior, ...
        innovation, ...
        measurementJacobian, ...
        measurementNoiseCovariance, ...
        gateThreshold);
    
    % Reject outlier without modifying state or covariance
    if ~accepted
        nominalStatePosterior = nominalStatePrior;
        errorCovariancePosterior = errorCovariancePrior;
        return;
    end
    
    nominalStatePosterior = state_injection(nominalStatePrior, errorCorrection);
    
    resetJacobian = compute_error_reset_jacobian(errorCorrection);
    
    errorCovariancePosterior = resetJacobian * errorCovarianceConditioned * resetJacobian';
    errorCovariancePosterior = 0.5 * (errorCovariancePosterior + errorCovariancePosterior');

end
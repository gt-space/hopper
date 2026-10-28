function [errorCorrection, errorCovariancePosterior, accepted] = ...
    kalman_measurement_update(errorCovariancePrior, innovation, measurementJacobian, measurementNoiseCovariance, gateThreshold)
% (Yuma's implementation | 10/23/26)
% Generic error-state Kalman measurement update
%
% Inputs:
% 
%   errorCovariancePrior        - Prior error covariance (19x19)
%   innovation                  - Predicted-minus-measured residual (mx1)
%   measurementJacobian         - Measurement error Jacobian (mx19)
%   measurementNoiseCovariance  - Measurement noise covariance (mxm)
%   gateThreshold               - Normalized innovation squared threshold
% 
% Outputs:
%
%   errorCorrection             - Estimated error correction state vector (19x1)
%   errorCovariancePosterior    - Measurement-conditioned (posterior) covariance (19x19)
%   accepted                    - Boolean flag if measurement passed innovation gating

    stateDim = size(errorCovariancePrior, 1);
    measDim = numel(innovation);

    innovation = innovation(:); % ensure column vector format

    % Validate dimensions
    assert(size(errorCovariancePrior, 2) == stateDim);
    assert(size(measurementJacobian, 1) == ...
        measDim && size(measurementJacobian,2) == stateDim);
    assert(all(size(measurementNoiseCovariance) == [measDim,measDim]));

    % Innovation covariance
    innovationCovariance = measurementJacobian * errorCovariancePrior * measurementJacobian' + measurementNoiseCovariance;
    innovationCovariance = 0.5 * (innovationCovariance + innovationCovariance'); % remove numerical asymmetry

    % REQUIRE positive-definite innovation covariance
    [choleskyFactor, flag] = chol(innovationCovariance, 'lower');
    if flag ~= 0
        warning('Innovation covariance must be positive definite!');
        errorCorrection = zeros(stateDim, 1);
        errorCovariancePosterior = errorCovariancePrior;
        accepted = false;
        return;
    end

    % Normalize Innovation Squared (NIS) for outlier rejection
    disp(size(innovation));
    disp(size(innovationCovariance));
    disp(size(choleskyFactor));
    disp(size(measurementJacobian));
    disp(size(measurementNoiseCovariance));
    whitenedInnovation = choleskyFactor \ innovation;
    nis = whitenedInnovation' * whitenedInnovation; % computes inner product (norm-squared)

    % Innovation gating (outlier rejection check)
    if nis > gateThreshold
        errorCorrection = zeros(stateDim, 1);
        errorCovariancePosterior = errorCovariancePrior;
        accepted = false;
        return;
    end

    accepted = true;

    % Kalman gain
    kalmanGain = (errorCovariancePrior * measurementJacobian') / innovationCovariance;

    % Estimated error correction vector
    errorCorrection = kalmanGain * innovation;
    % NOTE: This formula is under our chosen estimation-minus-truth error
    % convention, given that each sensor's H maps those error coordinates
    % into predicted-minus-measured innovation. That assumption remains to
    % be verified for all sensor models.

    
% Joseph-form covariance update
    identityMinusKH = eye(stateDim) - (kalmanGain * measurementJacobian);
    errorCovariancePosterior = identityMinusKH * errorCovariancePrior * identityMinusKH' + ...
        (kalmanGain * measurementNoiseCovariance * kalmanGain');

    errorCovariancePosterior = 0.5 * (errorCovariancePosterior + errorCovariancePosterior'); % remove numerical asymmetry

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update.m
% 
% Generic function to update a state estimate using a given sensor model
% and measurement.  
% 
% Author: David Meisinger
% Date: 10/18/26
% 
% Inputs:
%   x: state vector
%   P: covariance
%   model: function handle for sensor model
%   meas: measurement struct for the sensor
%   params: Any params needed for the sensor model
% 
% Outputs:
%   x_c: corrected state vector
%   P_c: corrected covariance
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
function [x_c, P_c] = update(x, P, model, meas, params)

    % Get the measurement innovation (estimate - meas) and Jacobian column
    [r, H] = model(x, meas, params);
    R = meas.R;
    S = H*P*H' + R;

    % chi-square innovation gate
    if r'*(S\r) > gate(numel(r))
        x_c = x; P_c = P;
        return;
    end
    
    K  = P*H'/S;
    dx = K*r;                              % 19x1 error estimate (est - truth)
    I_KH = eye(19) - K*H;
    P_c = I_KH*P*I_KH' + K*R*K';            % Joseph form
    P_c = 0.5*(P_c+P_c');
    x_c = inject(x, dx);                     % corrections subtract, since est - truth
end

%}
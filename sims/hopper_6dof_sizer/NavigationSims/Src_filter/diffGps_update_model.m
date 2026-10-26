% diffGps_update_model.m
%
% Function to compute the measurement innovation and Jacobian for the
% filter correct step. This function handle is passed into update.m for
% updating the state estimate.
%
% The measurement is the antenna 1 -> antenna 2 baseline in inertial NED:
%   h(x) = C_IB(q) * (r2 - r1)
% The baseline does not depend on position, so H is nonzero only in the
% attitude error block (10:12). Rotation about the baseline axis itself is
% unobservable from this measurement.
%
% Inputs:
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - differential GPS measurement struct
%   params              - Sensor parameters (gpsFrameBus params here)
%
% Outputs:
%   innovation          - Predicted-minus-measured residual (3x1)
%   measurementJacobian - Error-state measurement Jacobian H (3x19)
%
% Conventions:
%   - Quaternion represents inertial-to-body attitude
%   - Antenna positions r1, r2 are in body coordinates
%   - Attitude error is inertial-aligned (see state_injection.m), so a
%     true baseline b_true = (I - [dtheta]x) * b_predicted. With
%     innovation = predicted - measured, this gives H = -[b_predicted]x
function [innovation, measurementJacobian] = ...
    diffGps_update_model(nominalState, measurement, params)

    quaternionEstimate = nominalState(10:13);

    % Body-to-inertial rotation
    C_IB = quaternion_to_rotation(quaternionEstimate)';

    % Antenna 1 -> antenna 2 baseline in body coordinates
    baselineBody = params.r2(:) - params.r1(:);

    % Predicted baseline in NED
    predictedMeasurement = C_IB * baselineBody;

    % Innovation
    innovation = predictedMeasurement - measurement.basisVec(:);

    % Error-state measurement Jacobian
    measurementJacobian = zeros(3,19);
    measurementJacobian(:,10:12) = -skew(predictedMeasurement);

end

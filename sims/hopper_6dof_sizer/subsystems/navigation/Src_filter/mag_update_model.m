% mag_update_model.m
%
% Function to compute the measurement innovation and Jacobian for the
% filter correct step. This function handle is passed into update.m for
% updating the state estimate.
%
% The measurement is the Earth's magnetic field seen in body axes:
%   h(x) = C_IB(q)' * B_E
% where B_E is the field in inertial NED. Like the differential GPS
% baseline, it depends only on attitude, so H is nonzero only in the
% attitude error block (10:12). Rotation about B_E itself is unobservable
% from this measurement; with B_E ~60 deg from vertical it fixes heading.
%
% Inputs:
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - magnetometer measurement struct (bfield_meas,
%                         body axes; bfield_E, the field in NED)
%   params              - gate threshold only
%
% Outputs:
%   innovation          - Predicted-minus-measured residual (3x1)
%   measurementJacobian - Error-state measurement Jacobian H (3x19)
%
% Conventions:
%   - quaternion_to_rotation gives C_IB, body-to-inertial
%   - Attitude error is inertial-aligned (see state_injection.m and
%     diffGps_update_model.m): C_IB_true = (I - [dtheta]x) * C_IB_predicted,
%     so the true body field is m_predicted - C_IB' * [B_E]x * dtheta.
%     With innovation = predicted - measured this gives H = C_IB' * [B_E]x
function [innovation, measurementJacobian] = ...
    mag_update_model(nominalState, measurement, params) %#ok<INUSD>

    quaternionEstimate = nominalState(10:13);

    % Body-to-inertial rotation
    C_IB = quaternion_to_rotation(quaternionEstimate);

    % Earth field in NED (magBus.bfield_E, from the field model)
    fieldNED = measurement.bfield_E(:);

    % Predicted field in body coordinates
    predictedMeasurement = C_IB' * fieldNED;

    % Innovation
    innovation = predictedMeasurement - measurement.bfield_meas(:);

    % Error-state measurement Jacobian
    measurementJacobian = zeros(3,19);
    measurementJacobian(:,10:12) = C_IB' * skew(fieldNED);

end

function [innovation, measurementJacobian] = gyro_update_model(nominalState, measurement, params)
% (Yuma's implementation | 10/25/26)
% Gyroscope forward measurement model for the GMEKF
%
% Inputs:
%
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - Gyroscope measurement struct
%   params              - Sensor parameters (unused here)
%
% Outputs:
%
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state measurement Jacobian H (3x19)
%
% Conventions:
%   - Angular velocity and gyro bias stored in body coordinates
%   - Angular velocity error is body-frame additive
%   - Gyro bias error is inertial-aligned
%   - Quaternion represents inertial-to-body attitude
    
    angularVelocityEstimate = nominalState(7:9);
    quaternionEstimate = nominalState(10:13);
    gyroBiasEstimate = nominalState(17:19);

    predictedMeasurement = angularVelocityEstimate + gyroBiasEstimate;

    innovation = predictedMeasurement - measurement.angvel(:);

    % Inertial-to-body rotation
    C_BI = quaternion_to_rotation(quaternionEstimate)';

    measurementJacobian = zeros(3,19);

    % Angular velocity error
    measurementJacobian(:,7:9) = eye(3);

    % Inertial-aligned gyro bias error
    measurementJacobian(:,16:18) = C_BI;


end
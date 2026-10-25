function [innovation, measurementJacobian] = gps_update_model(nominalState, measurement, params)
% (Yuma's implementation | 10/25/26)
% GPS position/velocity forward measurement model for the GMEKF
%
% Inputs:
%
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - GPS measurement struct
%   params              - GPS parameters, including r_b (3x1)
%
% Outputs:
%
%   innovation          - Predicted-minus-measured residual (6x1)
%   measurementJacobian - Error-state measurement Jacobian H (6x19)
%
% Conventions:
%   - Position and velocity are in inertial NED
%   - Angular velocity is in body coordinates
%   - Quaternion represents inertial-to-body attitude
%   - Antenna lever arm r_b is in body coordinates


    positionEstimate = nominalState(1:3);
    velocityEstimate = nominalState(4:6);
    angularVelocityEstimate = nominalState(7:9);
    quaternionEstimate = nominalState(10:13);
    
    % Antenna lever arm
    r_B = params.r_b(:);

    % Body-to-inertial rotation
    C_IB = quaternion_to_rotation(quaternionEstimate)';

    % Antenna lever arm expressed in inertial coordinates
    leverArmInertial = C_IB * r_B;

    % Rotational contribution to antenna velocity
    rotationalVelocityInertial = C_IB * cross(angularVelocityEstimate, r_B);

    % Predicted GPS position and velocity
    predictedPosition = positionEstimate + leverArmInertial;
    predictedVelocity = velocityEstimate + rotationalVelocityInertial;
    predictedMeasurement = [predictedPosition; predictedVelocity];

    % Actual GPS measurement
    measuredGPS = [measurement.pos(:); measurement.vel(:)];

    % Innovation
    innovation = predictedMeasurement - measuredGPS;

    % Error-state measurement Jacobian
    measurementJacobian = zeros(6,19);

    % GPS position
    measurementJacobian(1:3,1:3) = eye(3);
    measurementJacobian(1:3,10:12) = -skew(leverArmInertial);

    % GPS velocity
    measurementJacobian(4:6,4:6) = eye(3);
    measurementJacobian(4:6,7:9) = -C_IB * skew(r_B);
    measurementJacobian(4:6,10:12) = -skew(rotationalVelocityInertial);

end
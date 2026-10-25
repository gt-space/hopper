<<<<<<< HEAD
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

=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gps_update_model.m
% 
% Computes the predicted GPS position measurement and its H matrix for the 
% filter.
% 
% Author: David Meisinger
% Date: 10/26/26
% 
% Inputs:
%   nominalState - Nominal state estimate (20x1)
%   measurement  - Barometer measurement struct
%   params       - Parameter struct. Contains r1 and r2 from gpsFrameBus
% 
% Outputs:
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state Jacobian (1x19)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [innovation, measurementJacobian] = gps_update_model(nominalState, measuremeant, params)

% Radii to antennae
r1 = params.r1;
r2 = params.r2;

% Convert quaternion to rotation matrix
R = quaternion_to_rotation(nominalState(10:13)); % Body to NED

% Predicted GPS position: position_estimated + R * r_b
positionPredicted = nominalState(1:3) + R * (r1+r2)/2;

% Predicted GPS velocity: velocity_estimated
velocityPredicted = nominalState(4:6);

% WIP: all code under this is TB's and needs adjustment. Description above
% is source of truth. Jacobian is suspect and needs to be expanded to
% include velocity.

% Measurement Jacobian H for GPS Position (3x19)
% State order: [delta_theta (3), nu (3), chi (3), beta (3), gamma (3), lambda (1)]
% For position measurement, chi affects it directly (identity), 
% and attitude error (delta_theta) affects the lever arm cross product term.
H_gps = zeros(3, 19);
% Derivative with respect to attitude error (delta_theta) via skew-symmetric matrix of R_m * r_b
skew_lever = [0, -(R(1,:)*r_b(3)), (R(1,:)*r_b(2)); ...
    (R(2,:)*r_b(3)), 0, -(R(2,:)*r_b(1)); ...
    -(R(3,:)*r_b(2)), (R(3,:)*r_b(1)), 0]; % approximation or exact skew
% Alternatively using cross-product skew:
% skew_lever = -skew_sym(R_m * r_b); 
H_gps(:, 1:3) = -skew_lever; % Attitude error coupling through lever arm
H_gps(:, 7:9) = eye(3);      % Position error (chi) coupling (identity)
>>>>>>> ef24f709a8dbddd6c79e22cd08112e04c7f92083
end
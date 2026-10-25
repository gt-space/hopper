% LIDAR forward measurement model for the GMEKF. Single beam model,
% params.lidarIndex selects which LIDAR to run the forward model for.
%
% Author: David Meisinger
% Date: 10/25/26
% 
% Inputs:
%   nominalState        - Nominal state estimate (20x1)
%   measurement         - LIDAR measurement struct
%   params              - Sensor parameters (lidarFrameBus in this case)
%
% Outputs:
%   innovation          - Predicted-minus-measured residual
%   measurementJacobian - Error-state measurement Jacobian H (1x19)
function [innovation, measurementJacobian] = lidar_update_model(nominalState, measurement, params)

% Load the LIDAR frame variables
idx = params.lidarIndex;
switch idx
    case 1, A = params.A1; r = params.r1;
    case 2, A = params.A2; r = params.r2;
    case 3, A = params.A3; r = params.r3;
    case 4, A = params.A4; r = params.r4;
    otherwise, error('lidar_update_model: lidarIndex must be 1-4.');
end

% Get the height and attitude of the vehicle
h = nominalState(3) - params.heightCOM;
q_BI = nominalState(10:13);

% Compute the rotation matrix from q
R = quaternion_to_rotation(q_BI);
n = R(:,3)';

% Compute the LIDAR distance estimate
D = n*A(:,1);
distanceEstimate = -(n*r + h)/D;

% Unusable geometry: beam (nearly) parallel to ground, or no hit ahead of sensor
if abs(D) < 1e-3 || d <= 0
    innovation                 = [];
    measurementJacobian        = zeros(0,19);
    return;
end

% Compute innovation
innovation = distanceEstimate - measurement.distance;

% Calculate measurement Jacobian using the GMEKF conventions
u = R*(r + distanceEstimate*A(:,1));                    % inertial COM -> ground-hit vector
measurementJacobian = zeros(1,19);
measurementJacobian(3)     = -1/D;                      % position z
measurementJacobian(10:12) = (cross([0 0 1]',u)/D)';    % attitude error (inertial frame)

end
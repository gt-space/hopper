%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gps_filter_model.m
% Computes the predicted GPS position measurement and its H matrix for the filter.
% 
% Inputs:
%   m_X: nominal position vector (3x1)
%   m_q: nominal attitude quaternion
%   gps_params: struct containing antenna lever arm (r_b)
% 
% Outputs:
%   h_gps: predicted GPS position measurement (3x1)
%   H_gps: measurement Jacobian matrix (3x19)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h_gps, H_gps] = gps_filter_model(m_X, m_q, gps_params)
r_b = gps_params.r_b;
% Convert quaternion to rotation matrix
R_m = quat2dcm(m_q)'; % Body to NED

% Predicted GPS position: p_true + C_bn^T * r_b[cite: 1]
h_gps = m_X + R_m * r_b;

% Measurement Jacobian H for GPS Position (3x19)
% State order: [delta_theta (3), nu (3), chi (3), beta (3), gamma (3), lambda (1)]
% For position measurement, chi affects it directly (identity), 
% and attitude error (delta_theta) affects the lever arm cross product term.
H_gps = zeros(3, 19);
% Derivative with respect to attitude error (delta_theta) via skew-symmetric matrix of R_m * r_b
skew_lever = [0, -(R_m(1,:)*r_b(3)), (R_m(1,:)*r_b(2)); ...
    (R_m(2,:)*r_b(3)), 0, -(R_m(2,:)*r_b(1)); ...
    -(R_m(3,:)*r_b(2)), (R_m(3,:)*r_b(1)), 0]; % approximation or exact skew
% Alternatively using cross-product skew:
% skew_lever = -skew_sym(R_m * r_b); 
H_gps(:, 1:3) = -skew_lever; % Attitude error coupling through lever arm
H_gps(:, 7:9) = eye(3);      % Position error (chi) coupling (identity)
end
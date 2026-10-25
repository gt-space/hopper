%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_update_model.m
% Computes the predicted LiDAR range and its Jacobian row for the filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho_hat, H_lidar] = lidar_update_model(m_X, m_q, lidar_params)
l = lidar_params.l;     % Lever arm (3x1)
d = lidar_params.d;     % Unit beam direction (3x1)
z_g = lidar_params.z_g; % Flat ground height
R_m = quat2dcm(m_q)';
e3 = [0; 0; 1];
D = e3' * R_m * d;
if abs(D) < 1e-6
    error('LiDAR beam is parallel to the ground or invalid intersection!');
end

% Predicted range equation[cite: 5]
numerator = e3' * (m_X + R_m * l) - z_g;
rho_hat = -numerator / D;

% Jacobian for MEKF/GMEKF as specified in your handoff hand-out[cite: 5]
% H_L^M = [-e3^T * R_m * [l + rho_hat*d]_x,  0_{1x3},  e3^T / D,  0_{1x3},  0_{1x3},  0_{1x1}]
H_lidar = zeros(1, 19);
vector_term = l + rho_hat * d;
% Skew symmetric cross product helper or manual expansion:
cross_skew = [0, -vector_term(3), vector_term(2); ...
    vector_term(3), 0, -vector_term(1); ...
    -vector_term(2), vector_term(1), 0];
H_lidar(1:3) = -e3' * R_m * cross_skew; % Attitude error block
H_lidar(7:9) = e3' / D;                 % Position error block
end
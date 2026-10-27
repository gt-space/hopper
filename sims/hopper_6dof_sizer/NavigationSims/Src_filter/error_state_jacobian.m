%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error_state_jacobian.m
% 
% Function to calculate the F and G matrices for the equation:
% xdot = F(t)x(t) + G(t)w(t) where w(t) is white noise.
% 
% Author: David Meisinger
% Date: 10/18/26
% 
% Inputs:
%   x: estimated state (20x1)
%   aux: struct with various params
%   J: moment of inertia tensor (kg*m^2)
% 
% Outputs:
%   F: State Jacobian
%   G: Process noise Jacobian
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, G] = error_state_jacobian(x, aux, J)

    % Init the outputs
    F = zeros(19);
    G = zeros(19,13);

    % De-package the state vector into useful variables
    angvel = x(7:9);

    % Get rotation matrix from attitude
    C = aux.A_IB05;

    % Get inertial acceleration
    a_I = aux.a_I;

    % Pull out the inertial angular velocity
    w_I = aux.w_I;

    % Calculate M
    M = J\(skew(J*angvel) - skew(angvel)*J);

    % Populate the F matrix
    F(1:3,   4:6)   = eye(3);                  % dp/dnu
    F(4:6,   10:12) = -skew(a_I);              % dnu/dalpha
    F(4:6,   13:15) = -eye(3);                 % dnu/dzeta
    % F(7:9,   7:9)   = skew(w_I) + C*M*C';      % deps/deps (Euler dynamics, in NED)
    F(7:9,   7:9)   = M;
    % F(10:12, 7:9)   = eye(3);                  % dalpha/deps
    F(10:12, 7:9)   = C;                  
    F(13:15, 13:15) = skew(w_I);               % dzeta/dzeta
    F(16:18, 16:18) = skew(w_I);               % dbeta/dbeta

    % Populate the G matrix

    G(4:6,   1:3)   = C;          % accel noise  -> nu
    G(7:9,   4:6)   = C/J;        % torque noise -> eps  (C*inv(J))
    G(7:9,4:6) = J \ eye(3);
    G(13:15, 7:9)   = C;          % accel-bias random walk
    G(16:18, 10:12) = C;          % gyro-bias random walk
    G(19,    13)    = 1;          % baro-bias random walk

end
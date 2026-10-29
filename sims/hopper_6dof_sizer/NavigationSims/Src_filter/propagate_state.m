%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagate_state.m
% 
% Function to take in the state x_n and IMU reading and propagate to x_n+1.
% Uses first-order Euler propagation (except second order for linear
% position) which may need to be upgraded later.
% 
% Author: David Meisinger
% Date: 10/14/26
% 
% Inputs:
%   xn: n state vector (20x1)
%   aux: struct with useful shared variables for the filter
%   accel_meas: IMU acceleration measured (3x1) (m/s^2) (body frame)
%   angvel_meas: IMU angular velocity measured (3x1) (deg/s) (body frame)
%   validity: if IMU measurement is valid
%   J: Inertia tensor, xyz aligning with body frame (kg*m^2)
%   g: acceleration of gravity (scalar)
%   dt: time delta to propagate (s)
% 
% Outputs:
%   xnp1: n+1 state vector (20x1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnp1 = propagate_state(xn, aux, accel_meas, angvel_meas, validity, J, g, dt)

    % First depackage the state vector into useful variables
    posn = xn(1:3);
    veln = xn(4:6);
    angveln = xn(7:9);
    qn = xn(10:13);
    abn = xn(14:16);
    gbn = xn(17:19);
    bbn = xn(20);

    J = J(1:3,1:3);

    % Predict state based on validity
    if validity
        % Acceleration of gravity vector in NED
        g_vec = [0; 0; g];
    
        % Update position and velocity
        accel_est = aux.A_IB*(accel_meas - abn) + g_vec;
        posnp1 = posn + veln*dt + accel_est/2*dt^2;
        velnp1 = veln + accel_est*dt;
        angvelnp1 = angvel_meas - abn;
    else
        % If invalid IMU measurement, dead-reckon based on previous state
        posnp1 = posn + veln*dt;
        velnp1 = veln;
        omegadot = -J \ (cross(angveln, J*angveln));
        angvelnp1 = angveln + omegadot*dt;
    end

    % Update quaternion using exponential with average angvel
    theta = (angveln+angvelnp1)/2*dt;
    ang = norm(theta);
    if ang < 10^-8
        dq = [1; 0.5*theta];
    else
        dq = [cos(ang/2); sin(ang/2)*theta/ang];
    end
    qnp1 = qmult(qn, dq);
    qnp1 = qnp1/norm(qnp1); % normalize the in case of numerical drift
    
    % Rebuild the propagated state vector
    xnp1 = [posnp1; velnp1; angvelnp1; qnp1; abn; gbn; bbn];
end
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
%   torque: control torque applied (3x1) (N*m) (body frame)
%   validity: if IMU measurement is valid
%   J: Inertia tensor, xyz aligning with body frame (kg*m^2)
%   g: acceleration of gravity (scalar)
%   dt: time delta to propagate (s)
% 
% Outputs:
%   xnp1: n+1 state vector (20x1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnp1 = propagate_state(xn, aux, accel_meas, torque, validity, J, g, dt)

    % First depackage the state vector into useful variables
    posn = xn(1:3);
    veln = xn(4:6);
    angveln = xn(7:9);
    qn = xn(10:13);
    abn = xn(14:16);
    gbn = xn(17:19);
    bbn = xn(20);

    % Predict state based on validity
    if validity
        % Acceleration of gravity vector in NED
        g_vec = [0; 0; g];
    
        % Update position and velocity
        accel_est = aux.A_IB*(accel_meas - abn) + g_vec;
        posnp1 = posn + veln*dt + accel_est/2*dt^2;
        velnp1 = veln + accel_est*dt;
    else
        % If invalid IMU measurement, dead-reckon based on previous state
        posnp1 = posn + veln*dt;
        velnp1 = veln;
    end

    % Update the rotational velocity using Heun (trapezoidal integration)
    f = @(om) J \ (torque - cross(om, J*om));
    k1 = f(angveln);
    angvelnp1 = angveln + 0.5*dt*(f(angveln) + f(angveln+dt*k1));
    angvelavg = 0.5*(angveln+angvelnp1);

    % Update quaternion using exponential with average angvel
    theta = angvelavg*dt;
    ang = norm(theta);
    if ang < 10^-8
        dq = [1; 0.5*theta];
    else
        dq = [cos(ang/2); sin(ang/2)*theta/ang];
    end
    dq_conj = [dq(1); -dq(2:4)];
    qnp1 = quatmult(dq_conj, qn);
    qnp1 = qnp1/norm(qnp1); % normalize the in case of numerical drift
    
    % Rebuild the propagated state vector
    xnp1 = [posnp1; velnp1; angvelnp1; qnp1; abn; gbn; bbn];
end
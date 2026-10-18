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
function xnp1 = propagate_state(xn, accel_meas, torque, validity, J, g, dt)

    % First depackage the state vector into useful variables
    posn = xn(1:3);
    veln = xn(4:6);
    angveln = xn(7:9);
    qn = xn(10:13);
    abn = xn(14:16);
    gbn = xn(17:19);
    bbn = xn(20);

    % Unpack quaternion
    w = qn(1); x = qn(2); y = qn(3); z = qn(4);

    % Predict state based on validity
    if validity
        % Acceleration of gravity vector in NED
        g_vec = [0; 0; g];
    
        % Update position and velocity
        A = [1-2*(y^2+z^2),  2*(x*y+w*z),    2*(x*z-w*y);
             2*(x*y-w*z),    1-2*(x^2+z^2),  2*(y*z+w*x);
             2*(x*z+w*y),    2*(y*z-w*x),    1-2*(x^2+y^2)];
        accel_est = A'*(accel_meas - abn) + g_vec;
        posnp1 = posn + veln*dt + accel_est/2*dt^2;
        velnp1 = veln + accel_est*dt;
    
        % Update the rotational velocity
        angveldot = inv(J)*(torque - cross(angveln, J*angveln));
        angvelnp1 = angveln + dt*angveldot;
    else
        % If invalid IMU measurement, dead-reckon based on previous state
        posnp1 = posn + veln*dt;
        velnp1 = veln;
        angvelnp1 = angveln;
    end

    % Update quaternion
    qthing = [-x -y -z;
               w -z  y;
               z  w -x;
              -y  x  w];
    qdot = 0.5*qthing*angvelnp1;
    qnp1 = qn + qdot*dt;
    qnp1 = qnp1/norm(qnp1); % normalize the quaternion
    
    % Rebuild the propagated state vector
    xnp1 = [posnp1; velnp1; angvelnp1; qnp1; abn; gbn; bbn];
end
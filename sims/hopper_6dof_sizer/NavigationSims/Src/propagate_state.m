%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagate_state.m
% 
% Function to take in the state x_n and IMU reading and propagate to x_n+1.
% Uses first-order Euler propagation (except second order for linear
% position) which may need to be upgraded later.
% 
% Author: David Meisinger
% Date: 10/11/26
% 
% Inputs:
%   xn: n state vector (20x1)
%   accel_meas: IMU acceleration measured (3x1) (m/s^2) (body frame)
%   angvel_meas: IMU angular vel measured (3x1) (radians/s) (body frame)
%   dt: time delta to propagate (s)
% 
% Outputs:
%   xnp1: n+1 state vector (20x1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnp1 = propagate_state(xn, accel_meas, angvel_meas, dt)
    
    % First depackage the state vector into useful variables
    posn = xn(1:3);
    veln = xn(4:6);
    qn = xn(10:13);
    abn = xn(14:16);
    gbn = xn(17:19);
    bbn = xn(20);

    % Acceleration of gravity vector in Atlanta
    g_vec = [0; 0; 9.795];

    % Update position and velocity
    w = qn(1); x = qn(2); y = qn(3); z = qn(4);
    A = [1-2*(y^2+z^2),  2*(x*y+w*z),    2*(x*z-w*y);
         2*(x*y-w*z),    1-2*(x^2+z^2),  2*(y*z+w*x);
         2*(x*z+w*y),    2*(y*z-w*x),    1-2*(x^2+y^2)];
    accel_est = A'*(accel_meas - abn) + g_vec;
    posnp1 = posn + veln*dt + accel_est/2*dt^2;
    velnp1 = veln + accel_est*dt;

    % Update the rotational state
    angvelnp1 = angvel_meas - gbn;
    qthing = [-x -y -z;
               w  z -y;
              -z  w  x;
               y -x  w];
    qdot = -0.5*qthing*angvelnp1;
    qnp1 = qn + qdot*dt;
    qnp1 = qnp1/norm(qnp1); % normalize the quaternion

    % Rebuild the propagated state vector
    xnp1 = [posnp1; velnp1; angvelnp1; qnp1; abn; gbn; bbn];

end
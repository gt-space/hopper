function xdot= hopper6dof2(t, x, T, delta_p, delta_y, F_rcs, m, g)
% 6-DOF hopper dynamics with time-varying mass & inertia
%% ==== STATES ====
X  = x(1); 
Y  = x(2); 
Z  = x(3);
u  = x(4);
v  = x(5);
w  = x(6);
P  = x(7); 
Q  = x(8);
R  = x(9);
q0 = x(10); 
q1 = x(11); 
q2 = x(12); 
q3 = x(13);

%% ==== CONSTANTS ====
m0  = 100;          % initial mass [kg]
mdot = 0.5;         % mass burn rate [kg/s]

D_m = 0.25;
r_m = D_m / 2;
L_m = 4;

%% ==== INERTIA (scaled with mass) ====
Ix0 = 0.5 * m0 * r_m^2;
Iy0 = (1/12) * m0 * (3*r_m^2 + L_m^2);
Iz0 = Iy0;

Ix = Ix0 * (m / m0);
Iy = Iy0 * (m / m0);
Iz = Iz0 * (m / m0);
Ixz = 0;

% Ix = 350;
% Iy = 350;
% Iz = 350;

%% ==== ROTATION MATRIX ====
Rbody2NED = [ ...
 q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
 2*(q1*q2+q0*q3),   q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
 2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   q0^2-q1^2-q2^2+q3^2 ];

%% ==== FORCES ====
Fx_T =  T * cosd(delta_p) * cosd(delta_y);
Fy_T = -T * cosd(delta_p) * sind(delta_y);
Fz_T = -T * sind(delta_p);

F_g_body = Rbody2NED.' * [0; 0; m*g];

Fx = Fx_T + F_g_body(1);
Fy = Fy_T + F_g_body(2);
Fz = Fz_T + F_g_body(3);

udot = Fx/m - Q*w + R*v;
vdot = Fy/m - R*u + P*w;
wdot = Fz/m - P*v + Q*u;

%% ==== MOMENTS ====
d = 1;   
r = [-d; 0; 0]; % In body frame x axes is along the centerline 
F_T = [Fx_T; Fy_T; Fz_T];
M_TVC = cross(r, F_T);

n_rcs = 2;
L = M_TVC(1) + r_m * F_rcs * n_rcs; % Roll moment 
M = M_TVC(2);                                     % Pitch moment 
N = M_TVC(3);                                       % Yaw moment 


% Rotational EOM
% p_dot
pdot = ( Iz*L + Ixz*N ...
       + (Ixz^2 - Iz*Iy)*Q*R ...
       + Ixz*(Ix - Iy)*P*Q ) ...
       / (Ix*Iz - Ixz^2);

% q_dot
qdot = ( M ...
       + (Iz - Ix)*P*R ...
       - Ixz*(P^2 - R^2) ) / Iy;

% r_dot
rdot = ( Ix*N + Ixz*L ...
       + (Ixz^2 - Ix*Iy)*P*Q ...
       + Ixz*(Iy - Iz)*Q*R ) ...
       / (Ix*Iz - Ixz^2);

%% ==== QUATERNION KINEMATICS ====
q0dot = -0.5*(q1*P + q2*Q + q3*R);
q1dot =  0.5*(q0*P + q2*R - q3*Q);
q2dot =  0.5*(q0*Q + q3*P - q1*R);
q3dot =  0.5*(q0*R + q1*Q - q2*P);

%% ==== POSITION KINEMATICS ====
pos_dot = Rbody2NED * [u; v; w];

%% ==== PACK ====
xdot = [
 pos_dot;
 udot; vdot; wdot;
 pdot; qdot; rdot;
 q0dot; q1dot; q2dot; q3dot ];

end
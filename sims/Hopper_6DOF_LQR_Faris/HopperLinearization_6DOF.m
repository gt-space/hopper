function   [A_num, B_num] = HopperLinearization_6DOF(x_val, u_val, m)
%% STATES & INPUTS
syms X Y Z u  v w P Q R q0 q1 q2 q3
% ========= STATES =============
x_linearization = [  X ; Y ; Z ;
                                u  ; v ; w ;
                                P  ; Q; R ;
                                q0 ; q1; q2 ; q3];

syms T delta_p delta_y F_rcs
u_linearization = [T; delta_p; delta_y; F_rcs];

%% Dynamics:
% ==== CONSTANT  ====
m0  = 100;          % initial mass [kg]
D_m = 0.25;
r_m = D_m/2;
L_m = 4;
g   = 9.81;

% ==== INERTIA (scaled with mass) ====
Ix0 = 0.5 * m0 * r_m^2;
Iy0 = (1/12) * m0 * (3*r_m^2 + L_m^2);
Iz0 = Iy0;

Ix = Ix0 * (m / m0);
Iy = Iy0 * (m / m0);
Iz = Iz0 * (m / m0);
Ixz = 0;

d = 1;   

% ==== ROTATION MATRIX ====
Rbody2NED = [ ...
 q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
 2*(q1*q2+q0*q3),   q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
 2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   q0^2-q1^2-q2^2+q3^2 ];

% ==== FORCES ====
Fx_T =  T * cosd(delta_p) * cosd(delta_y);
Fy_T = -T * cosd(delta_p) * sind(delta_y);
Fz_T = -T * sind(delta_p);

F_g_body = Rbody2NED.' * [0; 0; m*g];

Fx = Fx_T + F_g_body(1);
Fy = Fy_T + F_g_body(2);
Fz = Fz_T + F_g_body(3);

% ==== VELOCITY ====

udot = Fx/m - Q*w + R*v;
vdot = Fy/m - R*u + P*w;
wdot = Fz/m - P*v + Q*u;

% ==== MOMENTS ====
r = [-d; 0; 0];
F_T = [Fx_T; Fy_T; Fz_T];
M_TVC = cross(r, F_T);

n_rcs = 2;
L = M_TVC(1) + r_m * F_rcs * n_rcs;
M = M_TVC(2);
N = M_TVC(3);

pdot = (Iz*L + (Iz - Iy)*Q*R) / Ix;
qdot = (M + (Iz - Ix)*P*R) / Iy;
rdot = (Ix*N + (Ix - Iy)*P*Q) / Iz;

% ==== QUATERNION KINEMATICS ====
q0dot = -0.5*(q1*P + q2*Q + q3*R);
q1dot =  0.5*(q0*P + q2*R - q3*Q);
q2dot =  0.5*(q0*Q + q3*P - q1*R);
q3dot =  0.5*(q0*R + q1*Q - q2*P);

% ==== POSITION KINEMATICS ====
pos_dot = Rbody2NED * [u; v; w];

% ==== PACK ====
f = [
 pos_dot;
 udot; vdot; wdot;
 pdot; qdot; rdot;
 q0dot; q1dot; q2dot; q3dot ];
%% JACOBIANS: A and B
A = jacobian(f, x_linearization);      
B = jacobian(f, u_linearization);  


A_fun = matlabFunction(A,"Vars",{x_linearization, u_linearization});
B_fun = matlabFunction(B,"Vars",{x_linearization, u_linearization});


A_num = A_fun(x_val, u_val);
B_num = B_fun(x_val, u_val);

end

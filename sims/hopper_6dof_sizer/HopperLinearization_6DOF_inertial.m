function [A,B,f0] = HopperLinearization_6DOF_inertial(x, u, params)
% x = [X Y Z Vx Vy Vz P Q R q0 q1 q2 q3]'
% u = [T delta_p delta_y F_rcs]'   (radians)
% params.m, params.I (3x3), params.g, params.d, params.r_m

f0 = hopper_f(x,u,params);

h = 1e-12;                     % complex-step size
nx = numel(x); nu = numel(u);
A  = zeros(nx,nx);
B  = zeros(nx,nu);  

for i = 1:nx
    dx = zeros(nx,1); dx(i) = 1;
    fi = hopper_f(x + 1i*h*dx, u, params);
    A(:,i) = imag(fi)/h;
end

for j = 1:nu
    du = zeros(nu,1); du(j) = 1;
    fj = hopper_f(x, u + 1i*h*du, params);
    B(:,j) = imag(fj)/h;
end
end

function f = hopper_f(x,u,p)
X=x(1); Y=x(2); Z=x(3);
Vx=x(4); Vy=x(5); Vz=x(6);
P=x(7); Q=x(8); R=x(9);
q = x(10:13);   % [q0 q1 q2 q3]'

T       = u(1);
delta_p = u(2);
delta_y = u(3);
Tor_rcs   = u(4);

m = p.m; g = p.g; I = p.I; d = p.d; 

Rbn = Rbody2NED(q);


%check this
% Thrust in body (radians)
Fx_T =  T * cos(delta_p) * cos(delta_y);
Fy_T = -T * cos(delta_p) * sin(delta_y);
Fz_T = -T * sin(delta_p);
F_Tb = [Fx_T; Fy_T; Fz_T];

% Inertial translational dynamics (NED gravity)
F_n  = Rbn * F_Tb + [0;0;m*g];
Vdot = (1/m) * F_n;
posdot = [Vx;Vy;Vz];

% Moments in body
r_tvc = [-d;0;0];
M_tvc = cross(r_tvc, F_Tb);
M_body = [(M_tvc(1) +Tor_rcs);
          M_tvc(2);
          M_tvc(3)];

omega = [P;Q;R];
omegadot = I \ (M_body - cross(omega, I*omega));

% Quaternion kinematics
qdot = 0.5 * Omega(omega) * q;

f = [posdot;
     Vdot;
     omegadot;
     qdot];
end

function R = Rbody2NED(q)
q0=q(1); q1=q(2); q2=q(3); q3=q(4);
R = [ ...
 q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2);
 2*(q1*q2+q0*q3),   q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);
 2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   q0^2-q1^2-q2^2+q3^2 ];
end

function Om = Omega(w)
P=w(1); Q=w(2); R=w(3);
Om = [ 0  -P  -Q  -R;
       P   0   R  -Q;
       Q  -R   0   P;
       R   Q  -P   0 ];
end
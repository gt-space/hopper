function [A,B,f0] = HopperLinearization_6DOF(x,u,params)

f0 = hopper_f(x,u,params);

h  = 1e-20;

nx = numel(x);
nu = numel(u);

A = zeros(nx,nx);
B = zeros(nx,nu);

for i = 1:nx

    dx = zeros(nx,1);
    dx(i) = 1;

    f = hopper_f(x + 1i*h*dx , u , params);

    A(:,i) = imag(f)/h;

end

for j = 1:nu

    du = zeros(nu,1);
    du(j) = 1;

    f = hopper_f(x , u + 1i*h*du , params);

    B(:,j) = imag(f)/h;

end

end


function f = hopper_f(x,u,p)

%% STATE

pos = x(1:3);
vel = x(4:6);
omega = x(7:9);
q = x(10:13);

%% INPUT

T       = u(1);
delta_p = u(2);
delta_y = u(3);
F_rcs   = u(4);

%% PARAMETERS

m = p.m;
g = p.g;
I = p.I;
d = p.d;
r_m = p.r_m;

%% ROTATION

Rbn = quat2rot(q);

%% THRUST VECTOR (BODY)

Fx =  T*cos(delta_p)*cos(delta_y);
Fy = -T*cos(delta_p)*sin(delta_y);
Fz = -T*sin(delta_p);

Ftb = [Fx;Fy;Fz];

%% TRANSLATION

Fn = Rbn*Ftb + [0;0;m*g];

Vdot = Fn/m;

posdot = vel;

%% MOMENTS

rtvc = [-d;0;0];

Mtvc = cross(rtvc,Ftb);

Mbody = [ ...
Mtvc(1) + r_m*F_rcs
Mtvc(2)
Mtvc(3)];

omegadot = I\(Mbody - cross(omega,I*omega));

%% QUATERNION

qdot = 0.5*Omega(omega)*q;

%% OUTPUT

f = [posdot;Vdot;omegadot;qdot];

end


function R = quat2rot(q)

q0=q(1); q1=q(2); q2=q(3); q3=q(4);

R = [ ...
 q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2)
 2*(q1*q2+q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3-q0*q1)
 2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) q0^2-q1^2-q2^2+q3^2];

end


function Om = Omega(w)

P=w(1); Q=w(2); R=w(3);

Om = [ ...
0 -P -Q -R
P 0 R -Q
Q -R 0 P
R Q -P 0];

end
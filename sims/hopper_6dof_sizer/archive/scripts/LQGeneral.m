%clear all;
close all; clc;

clear A; clear B;


% x y z vx vy vz P Q R q0 q1 q2 q3
Q =  [1e3 0 0 0 0 0 0 0 0 0 0 0 0; % x
    0 1e3 0 0 0 0 0 0 0 0 0 0 0; % y
    0 0 1e13 0 0 0 0 0 0 0 0 0 0; % z
    0 0 0 1e3 0 0 0 0 0 0 0 0 0; % vx
    0 0 0 0 1e3 0 0 0 0 0 0 0 0; % vy
    0 0 0 0 0 5e9 0 0 0 0 0 0 0; % vz
    0 0 0 0 0 0 1e6 0 0 0 0 0 0; % P
    0 0 0 0 0 0 0 1e6 0 0 0 0 0; % Q
    0 0 0 0 0 0 0 0 1e6 0 0 0 0; % R
    0 0 0 0 0 0 0 0 0 100 0 0 0; % q0
    0 0 0 0 0 0 0 0 0 0 100 0 0; % q1
    0 0 0 0 0 0 0 0 0 0 0 100 0; % q2
    0 0 0 0 0 0 0 0 0 0 0 0 100]./10000; % q3

% T delta rcs
R = [3e8 0 0 0;
   0 8e5 0 0;
   0 0 8e5 0
   0 0 0 10000]./10000;

params = struct();
% constants
params.g = 9.80665;     % m/s^2  (Down is +Z in NED)
% params.m = linspace(OUT.Vehicle.WetMass,100, length(t));         % kg     (constant mass for now)
% geometry / actuator
params.d  = cg_init - engine_cg;       % m   thrust application x-offset in BODY frame (r = [-d;0;0])
params.r_m = 0.125;     % m   moment arm for RCS term (your r_m)
params.n_rcs = 1;       % unitless (if you keep it)


% inertia (constant, BODY frame about CG)
Ix = 1;   %[kg*m^2]
Iy = 41; 
Iz = 41; 
Ixy = 0; Ixz = 0; Iyz = 0;
params.I = [ Ix  -Ixy -Ixz;
           -Ixy  Iy  -Iyz;
           -Ixz -Iyz  Iz ];
% Preallocate (optional but nice)
nx = 13;
nu = 4;



Target_Trajectory = [0 0 0 0 0 0 0 0 0 sqrt(2)/2 0 sqrt(2)/2 0]';



% Linearize along the full target trajectory

   x = Target_Trajectory;   % 13x1
   u = [114*9.81; .5; .5; .5 ];             % 4x1
   [A,B] = HopperLinearization_6DOF_LQR(x, u, params);
 %  Fref(i,:) = f_i.';

Co = ctrb(A,B)
unco = length(A) - rank(Co)

N = zeros(13,4);
[K,S,P] = lqr(A,B,Q,R,N);



%%

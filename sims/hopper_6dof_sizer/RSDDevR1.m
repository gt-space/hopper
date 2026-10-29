% State order:
% [X Y Z Vx Vy Vz P Q R p1 p2 p3]'
%
% Input order:
% [T delta_p delta_y Tor_rcs]'
%
% Integral error convention:
% xi_dot = position_command - position

Qdiags  = [0 0 0 10 10 10 10 10 10 1e6 1e6 1e6];
QIdiags = [10 10 500];

Q  = diag(Qdiags);
QI = diag([Qdiags QIdiags]);

Rdiags = [1 1e8 1e8 1e8];
R = diag(Rdiags);

%% Parameters

params = struct();
params.g = 9.81;
params.d = cg_init - engine_cg;
params.r_m = 0.125;
params.n_rcs = 1;

% MoI_init must be the ORIGINAL inertia in the OLD body axes.
% Do not apply this transformation twice.
S = [ 0  0  1;
    0  1  0;
    -1  0  0];

params.I = S * MoI_init * S.';

% Use this same transformed inertia in the Simulink 6DOF block.
I_sim = params.I;

m_profile = linspace(OUT.Vehicle.WetMass, ...
    OUT.Vehicle.DryMass, length(t));
m = m_profile(1);

%% Hover equilibrium

p = zeros(3,1);

xbar = [zeros(9,1); p];
ubar = [m*params.g; 0; 0; 0];

% Preserve your original variable names.
x = xbar;
u = ubar;

[A,B] = HopperLinearization_lqi(xbar,ubar,params,m);

Cpos = [eye(3), zeros(3,9)];

%idek bro
A_i = [ 0 1;
    0 0 ];

B_i = 

A_a = [ A,     zeros(12,1);
          -B_i*C,  A_i];

B_a = [B; zeros(3,4)];




CrossTerm = zeros(15,4);

% RAW MATLAB gains:
% u = ubar - Klqi*[x-xbar; xi]
[Klqi,~,~] = lqr(A_a,B_a,QI,R,CrossTerm)

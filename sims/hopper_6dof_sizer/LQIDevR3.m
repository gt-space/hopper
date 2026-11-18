% State order:
% [X Y Z Vx Vy Vz P Q R p1 p2 p3]'

% Input order:
% [T delta_p delta_y Tor_rcs]'

% Integral error convention:
% xi_dot = position_command - position

Qdiags  = [1 1 1 10 10 100 10 10 10 1e5 1e5 1e5];
QIdiags = [100 100 1000000000];
QRamp = [1000000];

Q = diag([Qdiags QIdiags QRamp]);

Rdiags = [100 1e8 1e8 1e5];
R = diag(Rdiags);

%% Parameters

params = struct();
params.g = 9.81;
params.d = cg_init - engine_cg;
params.r_m = 0.125;
params.n_rcs = 1;

%change to new body axis
S = [ 0  0  1;
    0  1  0;
    -1  0  0];

params.I = S * MoI_init * S.';

I_sim = params.I;

m_profile = linspace(OUT.Vehicle.WetMass, ...
    OUT.Vehicle.DryMass, length(t));
m = m_profile(1);

%% Hover equilibrium

p = zeros(3,1);

x0 = [zeros(9,1); p];
u0 = [m*params.g; 0; 0; 0];

[A,B] = HopperLinearization_lqi(x0,u0,params,m);

%% Integral augmentation: xi_dot = r - y

C = [eye(3), zeros(3,9)];
D = zeros(3,4);

% Internal states:
% xi = [etaX; etaY; etaZ; rhoZ]

A_i = zeros(4);
A_i(4,3) = 1;       % rhoZ_dot = etaZ

B_i = [eye(3);
       zeros(1,3)];

n_i = 4; 

A_a = [A,       zeros(12,n_i);
       -B_i*C,  A_i];

B_a = [B;
       -B_i*D];
                 
CrossTerm = zeros(16,4);

[Klqi,~,~] = lqr(A_a,B_a,Q,R,CrossTerm)



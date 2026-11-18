% State order:
% [X Y Z Vx Vy Vz P Q R p1 p2 p3]'

% Input order:
% [T delta_p delta_y Tor_rcs]'

% Integral error convention:
% xi_dot = position_command - position

Qdiags  = [0 0 0 10 10 100 10 10 10 1e6 1e6 1e6];
QIdiags = [10 10 100000000];
QOscdiags = [5000 5000];

Q = diag([Qdiags QIdiags QOscdiags]);

Rdiags = [1 1e7 1e7 1e8];
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

% Three tracked outputs: X, Y, Z
Cpos = [eye(3), zeros(3,9)];
C = Cpos;
D = zeros(3,4);

% Angular frequency [rad/s].
omega_1 = 2*pi/54;

% Internal states:
% xi = [integral(eX); integral(eY); integral(eZ); z1; z2]

A_iPrime = [0,       -omega_1;
            omega_1, 0];

A_i = blkdiag(zeros(3), A_iPrime);

% e = [eX; eY; eZ]
% First three states integrate their respective errors.
% Last two states receive only eZ.
B_i = [eye(3);
       0 0 1;
       0 0 1];

n_i = 5;

A_a = [A,       zeros(12,n_i);
       -B_i*C,  A_i];

B_a = [B;
       -B_i*D];

CrossTerm = zeros(17,4);

[Klqi,~,poles] = lqr(A_a,B_a,Q,R,CrossTerm)



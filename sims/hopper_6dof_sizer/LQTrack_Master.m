%clear all;
close all;

load('thrust.mat');
APMAT = load("AP_T.mat");
T_profileEXCEL = readtable("ThrustProfileAP.xlsx");
T_profileEXCEL = table2array(T_profileEXCEL);

Trajectory2

%% -------------------------
% REDUCED 12-STATE WEIGHTS
% states = [x y z vx vy vz P Q R dth1 dth2 dth3]
% -------------------------

Sf = [1e3 0 0 0 0 0 0 0 0 0 0 0;   % x
      0 1e3 0 0 0 0 0 0 0 0 0 0;   % y
      0 0 1e10 0 0 0 0 0 0 0 0 0;  % z
      0 0 0 1e3 0 0 0 0 0 0 0 0;   % vx
      0 0 0 0 1e3 0 0 0 0 0 0 0;   % vy
      0 0 0 0 0 1e3 0 0 0 0 0 0;   % vz
      0 0 0 0 0 0 1e5 0 0 0 0 0;   % P
      0 0 0 0 0 0 0 1e5 0 0 0 0;   % Q
      0 0 0 0 0 0 0 0 1e5 0 0 0;   % R
      0 0 0 0 0 0 0 0 0 1e6 0 0;   % dth1
      0 0 0 0 0 0 0 0 0 0 1e5 0;   % dth2
      0 0 0 0 0 0 0 0 0 0 0 1e5]./10000;  % dth3

Q =  [1e3 0 0 0 0 0 0 0 0 0 0 0;      % x
      0 1e3 0 0 0 0 0 0 0 0 0 0;      % y
      0 0 1.7e13 0 0 0 0 0 0 0 0 0;   % z
      0 0 0 1e3 0 0 0 0 0 0 0 0;      % vx
      0 0 0 0 1e3 0 0 0 0 0 0 0;      % vy
      0 0 0 0 0 5e9 0 0 0 0 0 0;      % vz
      0 0 0 0 0 0 1e5 0 0 0 0 0;      % P
      0 0 0 0 0 0 0 1e5 0 0 0 0;      % Q
      0 0 0 0 0 0 0 0 1e5 0 0 0;      % R
      0 0 0 0 0 0 0 0 0 1e6 0 0;      % dth1
      0 0 0 0 0 0 0 0 0 0 1e4 0;      % dth2
      0 0 0 0 0 0 0 0 0 0 0 1e4]./10000;   % dth3

% T delta rcs
R = [3e7 0 0 0;
     0 8e6 0 0;
     0 0 8e6 0;
     0 0 0 1e6]./10000;

t = z_ts.Time(:,1);

AP = APMAT.ans;
t_old = double(AP.Time);
y_old = squeeze(AP.Data);
[t_old_unique, idx] = unique(t_old, 'stable');
y_old_unique = y_old(idx);

T_profile = T_ts.Data(:,1);

delta_y_profile = zeros(length(t),1);
delta_z_profile = zeros(length(t),1);
rcsp = zeros(length(t),1);

unom = [T_profile, delta_y_profile, delta_z_profile, rcsp];  % Nx4

tf = t(end);
t0 = 0;

x_trajectory = zeros(1,length(t));
y_trajectory = zeros(1,length(t));
z_trajectory = -z_ts.Data(:,1)';
z_trajectory(end) = 0;

vx_trajectory = zeros(1,length(t));
vy_trajectory = zeros(1,length(t));
vz_trajectory = -vz_ts.Data(:,1)';
vz_trajectory(end) = 0;
vz_trajectory(end-1) = 0;

P_trajectory = zeros(1,length(t));
Q_trajectory = zeros(1,length(t));
R_trajectory = zeros(1,length(t));

% full quaternion nominal trajectory still needed for linearization point
q0_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q1_trajectory = zeros(1,length(t));
q2_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q3_trajectory = zeros(1,length(t));

% ---- full 13-state nominal trajectory for linearization only
Target_Trajectory_full = [x_trajectory; y_trajectory; z_trajectory; ...
                          vx_trajectory; vy_trajectory; vz_trajectory; ...
                          P_trajectory; Q_trajectory; R_trajectory; ...
                          q0_trajectory; q1_trajectory; q2_trajectory; q3_trajectory ]';

% ---- reduced 12-state trajectory for gain solver
% last 3 states are reduced attitude-error coordinates, so nominal = 0
dth1_traj = zeros(1,length(t));
dth2_traj = zeros(1,length(t));
dth3_traj = zeros(1,length(t));

Target_Trajectory = [x_trajectory; y_trajectory; z_trajectory; ...
                     vx_trajectory; vy_trajectory; vz_trajectory; ...
                     P_trajectory; Q_trajectory; R_trajectory; ...
                     dth1_traj; dth2_traj; dth3_traj ]';

params = struct();
params.g = 9.80665;

params.d  = cg_init - engine_cg;
params.r_m = 0.125;
params.n_rcs = 1;

m_profile = linspace(OUT.Vehicle.WetMass,100,length(t));

Ix = 1;
Iy = 41;
Iz = 41;
Ixy = 0; Ixz = 0; Iyz = 0;

params.I = [ Ix  -Ixy -Ixz;
            -Ixy  Iy  -Iyz;
            -Ixz -Iyz  Iz ];

%% -------------------------
% REDUCED MODEL SIZES
% -------------------------
nx = 12;
nu = 4;
N  = length(t);

A = zeros(nx,nx,N);
B = zeros(nx,nu,N);
Fref = zeros(N,nx);

% Linearize about full nominal trajectory, but store reduced model
for i = 1:N
    x_i = Target_Trajectory_full(i,:).';   % full 13x1 nominal state
    u_i = unom(i,:).';                     % 4x1

    [A(:,:,i), B(:,:,i), f_i] = HopperLinearization_errorq(x_i, u_i, params, m_profile(i));
    Fref(i,:) = f_i.';
end

%% target trajectory (keep same length as t)
tgrid = t;

GainTable = TrajectoryFollowingGains3(Sf,Q,R,Target_Trajectory,A,B,tgrid,unom,params,Fref);

%% gains
K1 = GainTable.K1;
K2grid = GainTable.K2';

alpha = 1.04;
tau = 3;

bias = alpha * unom(:,1) .* exp(-tgrid/tau);

K2grid_mod = K2grid;
K2grid_mod(:,1) = K2grid(:,1) - bias;

%% -------------------------
% K1 flatten size changes: 4 x 12 = 48
% -------------------------
K1flat = zeros(length(tgrid),48);

for i = 1:length(tgrid)
    K1flat(i,:) = reshape(K1(:,:,i).',1,48);
end

z_ref_ts  = timeseries(z_trajectory, t);
T_ref_ts  = timeseries(T_profile, t);
vz_ref_ts = timeseries(vz_trajectory, t);

figure();
hold on;
for i = 1:12
    plot(t, squeeze(K1(2,i,:)), 'DisplayName', sprintf('K1(1,%d)', i));
end
xlabel('t');
ylabel('K1 second row values');
title('K1(1,:) vs time');
legend('show');
grid on;

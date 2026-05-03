%clear all;
close all;


Trajectory2


% x y z vx vy vz P Q R q0 q1 q2 q3
Sf = [1e6 0 0 0 0 0 0 0 0 0 0 0 0; % x
    0 1e6 0 0 0 0 0 0 0 0 0 0 0; % y
    0 0 1e10 0 0 0 0 0 0 0 0 0 0; % z
    0 0 0 1e3 0 0 0 0 0 0 0 0 0; % vx
    0 0 0 0 1e3 0 0 0 0 0 0 0 0; % vy
    0 0 0 0 0 1e3 0 0 0 0 0 0 0; % vz
    0 0 0 0 0 0 1e6 0 0 0 0 0 0; % P
    0 0 0 0 0 0 0 1e6 0 0 0 0 0; % Q
    0 0 0 0 0 0 0 0 1e6 0 0 0 0; % R
    0 0 0 0 0 0 0 0 0 1e5 0 0 0; % q0
    0 0 0 0 0 0 0 0 0 0 1e5 0 0; % q1
    0 0 0 0 0 0 0 0 0 0 0 1e5 0; % q2
    0 0 0 0 0 0 0 0 0 0 0 0 1e5]./10000; % q3

Q =  [1e3 0 0 0 0 0 0 0 0 0 0 0 0; % x
    0 1e3 0 0 0 0 0 0 0 0 0 0 0; % y
    0 0 1.7e13 0 0 0 0 0 0 0 0 0 0; % z
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
R = [3e7 0 0 0;
   0 8e5 0 0;
   0 0 8e5 0
   0 0 0 10000]./10000;



t= z_ts.Time(:,1);

T_profile = T_ts.Data(:,1);

delta_y_profile = zeros(length(t),1);        % Nx1
delta_z_profile = zeros(length(t),1);        % Nx1
rcsp = zeros(length(t),1);                   % Nx1

Tprof2 = ones(length(t),1);
unom = [T_profile, delta_y_profile, delta_z_profile, rcsp];  % Nx4

tf = t(end);
t0=0;

x_trajectory = zeros(1,length(t));
y_trajectory = zeros(1,length(t));
z_trajectory = -z_ts.Data(:,1)';

vx_trajectory = zeros(1,length(t));
vy_trajectory = zeros(1,length(t));
vz_trajectory = diff(z_trajectory)./diff(t)';
vz_trajectory(end+1)=vz_trajectory(end);

P_trajectory = zeros(1,length(t));
Q_trajectory = zeros(1,length(t));
R_trajectory = zeros(1,length(t));

q0_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q1_trajectory = zeros(1,length(t));
q2_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q3_trajectory = zeros(1,length(t));


Target_Trajectory = [x_trajectory;y_trajectory;z_trajectory;...
   vx_trajectory;vy_trajectory; vz_trajectory;...
   P_trajectory; Q_trajectory; R_trajectory;
   q0_trajectory; q1_trajectory; q2_trajectory;q3_trajectory ]';

params = struct();
% constants
params.g = 9.80665;     % m/s^2  (Down is +Z in NED)
% params.m = linspace(OUT.Vehicle.WetMass,100, length(t));         % kg     (constant mass for now)
% geometry / actuator
params.d  = cg_init - engine_cg;       % m   thrust application x-offset in BODY frame (r = [-d;0;0])
params.r_m = 0.125;     % m   moment arm for RCS term (your r_m)
params.n_rcs = 1;       % unitless (if you keep it)
m_profile = linspace(OUT.Vehicle.WetMass,OUT.Vehicle.DryMass, length(t));

% inertia (constant, BODY frame about CG)
params.I = MoI_init;

% Preallocate (optional but nice)
nx = 13;
nu = 4;
N  = length(t);
A  = zeros(nx,nx,N);
B  = zeros(nx,nu,N);

% Linearize along the full target trajectory
for i = 1:N
   x_i = Target_Trajectory(i,:).';   % 13x1
   u_i = unom(i,:).';               % 4x1
   [A(:,:,i), B(:,:,i),f_i] = HopperLinearization_6DOF_inertial(x_i, u_i, params, m_profile(i));
   Fref(i,:) = f_i.';

end
%% target trajectory (keep same length as t)
% Use the cleaned mission time grid directly
tgrid = t;
% Call gains using consistent (unique) grid + consistent trajectory
%GainTable = TrajectoryFollowingGains(Sf,Q,R,Target_Trajectory,A,B,tgrid,T_profile);
GainTable = TrajectoryFollowingGains2(Sf,Q,R,Target_Trajectory,A,B,tgrid,unom,params,Fref);

K1 = GainTable.K1;
K2grid = GainTable.K2';   

alpha = 1.04;
tau = 3;   % seconds

bias = alpha * unom(:,1) .* exp(-tgrid/tau);

K2grid_mod = K2grid;
K2grid_mod(:,1) = K2grid(:,1) - bias;


K1flat = zeros(length(tgrid),52);

for i = 1:length(tgrid)
   K1flat(i,:) = reshape(K1(:,:,i).',1,52);
end

z_ref_ts = timeseries(z_trajectory', t);
T_ref_ts = timeseries(T_profile, t);
vz_ref_ts = timeseries(vz_trajectory', t);


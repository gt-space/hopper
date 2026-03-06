clear all; close all; clc;

load('thrust.mat');
load('vz.mat');
load('z.mat');

% x y z vx vy vz P Q R q0 q1 q2 q3

Sf = [1 0 0 0 0 0 0 0 0 0 0 0 0; % x
      0 1 0 0 0 0 0 0 0 0 0 0 0; % y
      0 0 100 0 0 0 0 0 0 0 0 0 0; % z
      0 0 0 1 0 0 0 0 0 0 0 0 0; % vx
      0 0 0 0 1 0 0 0 0 0 0 0 0; % vy
      0 0 0 0 0 100000000000 0 0 0 0 0 0 0; % vz
      0 0 0 0 0 0 1 0 0 0 0 0 0; % P
      0 0 0 0 0 0 0 1 0 0 0 0 0; % Q
      0 0 0 0 0 0 0 0 1 0 0 0 0; % R
      0 0 0 0 0 0 0 0 0 1 0 0 0; % q0
      0 0 0 0 0 0 0 0 0 0 1 0 0; % q1
      0 0 0 0 0 0 0 0 0 0 0 1 0; % q2
      0 0 0 0 0 0 0 0 0 0 0 0 1]; % q3

Q =  [1 0 0 0 0 0 0 0 0 0 0 0 0; % x
      0 1 0 0 0 0 0 0 0 0 0 0 0; % y
      0 0 100000000 0 0 0 0 0 0 0 0 0 0; % z
      0 0 0 1 0 0 0 0 0 0 0 0 0; % vx
      0 0 0 0 1 0 0 0 0 0 0 0 0; % vy
      0 0 0 0 0 10 0 0 0 0 0 0 0; % vz
      0 0 0 0 0 0 1 0 0 0 0 0 0; % P
      0 0 0 0 0 0 0 1 0 0 0 0 0; % Q
      0 0 0 0 0 0 0 0 1 0 0 0 0; % R
      0 0 0 0 0 0 0 0 0 1000 0 0 0; % q0
      0 0 0 0 0 0 0 0 0 0 1000 0 0; % q1
      0 0 0 0 0 0 0 0 0 0 0 1000 0; % q2
      0 0 0 0 0 0 0 0 0 0 0 0 1000]; % q3

% T delta
R = [10 0 0 0;
     0 1 0 0; 
     0 0 1 0 
     0 0 0 1];

t = thrust_mat(:,1);

% T_profile = thrust_mat(:,2)';
% 
% %delta_p and delta z mixup?
% delta_y_profile = zeros(1,length(t));
% delta_z_profile = zeros(1,length(t));
% rcsp = zeros(1,length(t));
% 
% unom = [T_profile,delta_y_profile,delta_z_profile,rcsp];

T_profile = thrust_mat(:,2);                 % Nx1  (column)
delta_y_profile = zeros(length(t),1);        % Nx1
delta_z_profile = zeros(length(t),1);        % Nx1
rcsp = zeros(length(t),1);                   % Nx1

unom = [T_profile, delta_y_profile, delta_z_profile, rcsp];  % Nx4

tf = t(end);
t0=0;


%theta = ones(1,length(t)).*pi/2;
%delta = zeros(1,length(t));

x_trajectory = zeros(1,length(t));
y_trajectory = zeros(1,length(t));
z_trajectory = -z_mat(:,2)';
vx_trajectory = zeros(1,length(t));
vy_trajectory = zeros(1,length(t));
vz_trajectory = -vz_mat(:,2)';
P_trajectory = zeros(1,length(t));
Q_trajectory = zeros(1,length(t));
R_trajectory = zeros(1,length(t));

%up?
q0_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q1_trajectory = zeros(1,length(t));
q2_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q3_trajectory = zeros(1,length(t));


%theta_trajectory = (pi/2).*ones(1,length(t));
%q_trajectory = zeros(1,length(t));

Target_Trajectory = [x_trajectory;y_trajectory;z_trajectory;...
    vx_trajectory;vy_trajectory; vz_trajectory;...
    P_trajectory; Q_trajectory; R_trajectory;
    q0_trajectory; q1_trajectory; q2_trajectory;q3_trajectory ]';

params = struct();

% constants
params.g = 9.80665;     % m/s^2  (Down is +Z in NED)
params.m = 114;         % kg     (constant mass for now)

% geometry / actuator
params.d   = 1.0;       % m   thrust application x-offset in BODY frame (r = [-d;0;0])
params.r_m = 0.125;     % m   moment arm for RCS term (your r_m)
params.n_rcs = 1;       % unitless (if you keep it)

% inertia (constant, BODY frame about CG)

Ix = 16.5;   %[kg*m^2]
Iy = 16.5;  
Iz = 16.5;  
Ixy = 0; Ixz = 0; Iyz = 0;

params.I = [ Ix  -Ixy -Ixz;
            -Ixy  Iy  -Iyz;
            -Ixz -Iyz  Iz ];

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
    [A(:,:,i), B(:,:,i)] = HopperLinearization_6DOF_inertial(x_i, u_i, params);
end


%% target trajectory (keep same length as t)

% Use the cleaned mission time grid directly
tgrid = t;

% Call gains using consistent (unique) grid + consistent trajectory
GainTable = TrajectoryFollowingGains(Sf,Q,R,Target_Trajectory,A,B,tgrid,T_profile);

%%

K1 = GainTable.K1;
K2grid = GainTable.K2';   % keep your Simulink convention (N×2)


K2grid = zeros(length(t),4);

%fix for 6dof
% Flatten K1 for Simulink (N×12)  

%4*13?
K1flat = zeros(length(tgrid),52);
for i = 1:length(tgrid)
    K1flat(i,:) = reshape(K1(:,:,i).',1,52);
end


z_ref_ts = timeseries(z_trajectory, t);
T_ref_ts = timeseries(T_profile, t);
vz_ref_ts = timeseries(vz_trajectory, t);
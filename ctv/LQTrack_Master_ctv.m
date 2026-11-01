%% LQTrack_Master_ctv
% NOT THE ACTIVE SETUP PATH. v2 is the model we fly, and sim_setup_ctv.m
% runs LinerizationMaster_ctv for it. This file is the v3 setup, kept for
% whenever v3 is picked up again. See the note at the end about what v3
% still needs beyond this.
%
% Builds the reference trajectory and the time-varying LQR gain table that
% hopper_6dof_NED_v3_ctv.slx actually consumes.
%
% WHY THIS FILE EXISTS
% The v3 model's Gain_Selector does
%     K1 = reshape(Kflat, 12, 4)'
% and reads a Constant block whose value is Target_Trajectory_full, i.e. it
% wants the REDUCED 12-state controller (4x12 gains -> 48-column K1flat) and
% a separate full 13-state reference for the attitude error computation.
%
% LinerizationMaster_ctv.m builds the FULL 13-state controller instead:
% Target_Trajectory (no _full) and a 52-column K1flat (4x13). That is what
% hopper_6dof_NED_v2_ctv.slx wants -- its Gain_Selector does
% reshape(Kflat, 13, 4). So LinerizationMaster_ctv.m is the v2 setup and is
% kept for that model; this file is the v3 setup.
%
% Ported from sims/hopper_6dof_sizer/LQTrack_Master.m with three changes:
%   - uses ctv's own cg/MoI (MoI_init, cg_init) instead of hardcoded inertia
%   - mass profile runs wet -> dry from OUT, instead of a hardcoded 100 kg
%   - (roll is kept as a 4th control input: CTV has roll control via four
%     small propeller motors, modelled by the EDF Dynamics subsystem)
%
% Produces, in the base workspace:
%   tgrid, Target_Trajectory_full (Nx13), Target_Trajectory (Nx12),
%   unom (Nx4), K1flat (Nx48), K2grid (Nx4), K2grid_mod (Nx4)

close all;

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

% T delta_p delta_y rcs
R = [3e7 0 0 0;
     0 8e6 0 0;
     0 0 8e6 0;
     0 0 0 1e6]./10000;

t = z_ts.Time(:,1);

T_profile = T_ts.Data(:,1);

delta_y_profile = zeros(length(t),1);
delta_z_profile = zeros(length(t),1);
rcsp            = zeros(length(t),1);

unom = [T_profile, delta_y_profile, delta_z_profile, rcsp];  % Nx4

tf = t(end);
t0 = 0;

x_trajectory = zeros(1,length(t));
y_trajectory = zeros(1,length(t));
z_trajectory = -z_ts.Data(:,1)';
z_trajectory(end) = 0;

vx_trajectory = zeros(1,length(t));
vy_trajectory = zeros(1,length(t));

% ctv's Trajectory2 does not publish vz_ts, so differentiate the position
% reference on its own grid (same approach as LinerizationMaster_ctv.m).
vz_trajectory = diff(z_trajectory)./diff(t)';
vz_trajectory(end+1)   = vz_trajectory(end);
vz_trajectory(end)     = 0;
vz_trajectory(end-1)   = 0;

P_trajectory = zeros(1,length(t));
Q_trajectory = zeros(1,length(t));
R_trajectory = zeros(1,length(t));

% full quaternion nominal trajectory still needed for linearization point
q0_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q1_trajectory = zeros(1,length(t));
q2_trajectory = ones(1,length(t)).*(sqrt(2)/2);
q3_trajectory = zeros(1,length(t));

% ---- full 13-state nominal trajectory (this is what the model reads)
Target_Trajectory_full = [x_trajectory; y_trajectory; z_trajectory; ...
                          vx_trajectory; vy_trajectory; vz_trajectory; ...
                          P_trajectory; Q_trajectory; R_trajectory; ...
                          q0_trajectory; q1_trajectory; q2_trajectory; q3_trajectory ]';

% ---- reduced 12-state trajectory for the gain solver
% last 3 states are reduced attitude-error coordinates, so nominal = 0
dth1_traj = zeros(1,length(t));
dth2_traj = zeros(1,length(t));
dth3_traj = zeros(1,length(t));

Target_Trajectory = [x_trajectory; y_trajectory; z_trajectory; ...
                     vx_trajectory; vy_trajectory; vz_trajectory; ...
                     P_trajectory; Q_trajectory; R_trajectory; ...
                     dth1_traj; dth2_traj; dth3_traj ]';

params = struct();
params.g   = 9.80665;
params.d   = cg_init - engine_cg;
params.r_m = 0.125;
params.n_rcs = 1;
params.I   = MoI_init;

m_profile = linspace(OUT.Vehicle.WetMass, OUT.Vehicle.DryMass, length(t));

%% -------------------------
% REDUCED MODEL SIZES
% -------------------------
nx = 12;
nu = 4;
N  = length(t);

A    = zeros(nx,nx,N);
B    = zeros(nx,nu,N);
Fref = zeros(N,nx);

% Linearize about the full nominal trajectory, but store the reduced model
for i = 1:N
    x_i = Target_Trajectory_full(i,:).';   % full 13x1 nominal state
    u_i = unom(i,:).';                     % 4x1

    [A(:,:,i), B(:,:,i), f_i] = ...
        HopperLinearization_errorq_ctv(x_i, u_i, params, m_profile(i));
    Fref(i,:) = f_i.';
end

%% target trajectory (keep same length as t)
tgrid = t;

GainTable = TrajectoryFollowingGains3(Sf,Q,R,Target_Trajectory,A,B,tgrid,unom,params,Fref);

K1 = GainTable.K1;
K2grid = GainTable.K2';

alpha = 1.04;
tau   = 3;

bias = alpha * unom(:,1) .* exp(-tgrid/tau);

K2grid_mod = K2grid;
K2grid_mod(:,1) = K2grid(:,1) - bias;

%% -------------------------
% K1 flatten size: 4 x 12 = 48, which is what the v3 Gain_Selector's
% reshape(Kflat, 12, 4)' expects.
% -------------------------
K1flat = zeros(length(tgrid),48);

for i = 1:length(tgrid)
    K1flat(i,:) = reshape(K1(:,:,i).',1,48);
end

z_ref_ts  = timeseries(z_trajectory, t);
T_ref_ts  = timeseries(T_profile, t);
vz_ref_ts = timeseries(vz_trajectory, t);

%% -------------------------
% STILL OUTSTANDING FOR v3 (not needed by v2):
% the v3 Slosh subsystem references IN.tanks.lateral_damping and
% IN.tanks.axial_damping, which ctv/inputs/mission_inputs.m leaves
% commented out and generateBaffleLUT_ctv does not assign. Uncomment those
% two lines (the sizer's copy has them, marked "MUST BE UNCOMMENTED")
% before trying to simulate v3.
% -------------------------

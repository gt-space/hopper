%clear all;
close all;


Trajectory2
ThreeSecTrajectory2


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
    0 0 0 0 0 0 0 0 0 1e4 0 0 0; % q0
    0 0 0 0 0 0 0 0 0 0 1e4 0 0; % q1
    0 0 0 0 0 0 0 0 0 0 0 1e4 0; % q2
    0 0 0 0 0 0 0 0 0 0 0 0 1e4]./10000; % q3
% Retuned 2026-10-20 on the rigid (slosh-off) vehicle: quaternion weights
% raised from 100 and TVC input weights (R below) raised 10x. The old gains
% drove the gimbal into its 15 deg/s rate limit and only flew with
% double-counted slosh inertia in the plant. These
% fly with the plant inertia scaled 0.8x-1.2x.

% T delta rcs
R = [3e7 0 0 0;
   0 8e6 0 0;
   0 0 8e6 0
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
params.d  = cg_init - IN.mount_cg;     % m   thrust line through the gimbal pivot, x-offset in BODY frame (r = [-d;0;0])
params.r_m = 0.125;     % m   moment arm for RCS term (your r_m)
params.n_rcs = 1;       % unitless (if you keep it)
% Propellant along the reference trajectory: ox/fuel flow from the
% reference thrust (mdot_lookup), so mass, CG and inertia follow the
% propellant actually burned rather than a straight line to dry mass.
% mdot_lookup.xlsx holds ox in column 2 and fuel in column 3
% (mdot_table_generation), but load_lookup names them the other way round
% and the sim's mdot Subsystem swaps its outports back: in the workspace,
% fu_mdot_data is the oxidizer flow and ox_mdot_data the fuel flow.
ox_mdot_ref = interp1(thrust_data, fu_mdot_data, T_profile, 'linear', 'extrap');
fu_mdot_ref = interp1(thrust_data, ox_mdot_data, T_profile, 'linear', 'extrap');
m_ox_ref = max(IN.propulsion.oxidizer_mass - cumtrapz(t, ox_mdot_ref), 0);
m_fu_ref = max(IN.propulsion.fuel_mass     - cumtrapz(t, fu_mdot_ref), 0);
m_profile = (OUT.Vehicle.DryMass + m_ox_ref + m_fu_ref).';

% inertia (constant, BODY frame about CG)
params.I = MoI_init;

% Preallocate (optional but nice)
nx = 13;
nu = 4;
N  = length(t);
A  = zeros(nx,nx,N);
B  = zeros(nx,nu,N);

% Schedule inertia and thrust moment arm with the propellant left along
% the reference (mass_properties_setup dry block + propellant columns);
% false -> constant full-tank values
if ~exist('lqr_schedule_inertia', 'var'), lqr_schedule_inertia = true; end

% Linearize along the full target trajectory
for i = 1:N
   x_i = Target_Trajectory(i,:).';   % 13x1
   u_i = unom(i,:).';               % 4x1
   params_i = params;
   if lqr_schedule_inertia
       [cg_i, params_i.I] = vehicle_mass_props(MP.DRY, MP.TG, m_ox_ref(i), m_fu_ref(i));
       params_i.d = cg_i(1) - IN.mount_cg;
   end
   [A(:,:,i), B(:,:,i),f_i] = HopperLinearization_6DOF_inertial(x_i, u_i, params_i, m_profile(i));
   Fref(i,:) = f_i.';

end
%% target trajectory (keep same length as t)
% Use the cleaned mission time grid directly
tgrid = t;
% Call gains using consistent (unique) grid + consistent trajectory
%GainTable = TrajectoryFollowingGains(Sf,Q,R,Target_Trajectory,A,B,tgrid,T_profile);
% The finite-horizon Riccati sweep starts from Sf at the final time. Over
% the last ~1.5 s that drove the TVC position gain up ~13x and the attitude
% gain to zero (K = R\B'*Sf only damps rates), so any drift at touchdown
% (e.g. from wind) was corrected with a large, unrestored tilt.
% Add a run-out to the reference: after the trajectory ends, keep
% descending at its final rate with hover thrust for lqr_runout_time
% seconds. The Riccati end effects then fall after touchdown, and a vehicle
% still above ground when the trajectory ends keeps descending instead of
% settling into a hover offset from a frozen final point. The run-out is
% part of the exported reference and gain table; the sim stops at ground
% contact, so the model StopTime must cover tgrid(end).
if ~exist('lqr_runout_time', 'var'), lqr_runout_time = 5; end   % s
dt_end = tgrid(end) - tgrid(end-1);
n_run  = round(lqr_runout_time / dt_end);
if n_run > 0
    tau    = dt_end * (1:n_run)';
    x_end  = Target_Trajectory(end,:);
    TT_run = repmat(x_end, n_run, 1);
    TT_run(:, 4:5) = 0;                              % no horizontal motion
    TT_run(:, 7:9) = 0;                              % no body rates
    TT_run(:, 6)   = x_end(6);                       % final descent rate (NED, down +)
    TT_run(:, 3)   = x_end(3) + x_end(6) * tau;      % keep descending
    u_run  = repmat([m_profile(end) * params.g, 0, 0, 0], n_run, 1);   % hover thrust
    [A_r, B_r, f_r] = HopperLinearization_6DOF_inertial(TT_run(1,:).', u_run(1,:).', params_i, m_profile(end));
    tgrid  = [tgrid; tgrid(end) + tau];
    Target_Trajectory = [Target_Trajectory; TT_run];
    unom   = [unom; u_run];
    A      = cat(3, A, repmat(A_r, 1, 1, n_run));
    B      = cat(3, B, repmat(B_r, 1, 1, n_run));
    Fref   = [Fref; repmat(f_r.', n_run, 1)];
end
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

% --- EXPORT EVERYTHING TO BASE WORKSPACE FOR SIMULINK ---
assignin('base', 'tgrid',                 tgrid);
assignin('base', 'unom',                  unom);
assignin('base', 'Target_Trajectory',     Target_Trajectory);
assignin('base', 'z_ref_ts',              z_ref_ts);
assignin('base', 'vz_ref_ts',             vz_ref_ts);
assignin('base', 'K1flat',                K1flat);
assignin('base', 'K2grid',                K2grid);
assignin('base', 'K2grid_mod',            K2grid_mod);

% Check if your model uses T_ts or T_ref_ts (pushing both covers all bases)
if exist('T_ts', 'var')
    assignin('base', 'T_ts',              T_ts);
end
if exist('T_ref_ts', 'var')
    assignin('base', 'T_ref_ts',          T_ref_ts);
end

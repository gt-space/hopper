close all

% Simulation Parameters
dt = 0.1; % Time step (s)
T = 5; % Time Span (s)
N = T/dt + 1; % Number of simulated states

% Vehicle Parameters
vehicle.g = 9.81; % Gravity (m/s^2)
vehicle.Tmax = 1556.88; % max thrust [N]
vehicle.Isp = 260; % seconds
vehicle.mr = 3.8; % O/F Mixture Ratio
vehicle.n2o_mass_0 = 40; % Initial Propellant Masses 
vehicle.ipa_mass_0 = 9; 
vehicle.n2o_tank_height = 0.5; % Propellant tank heights
vehicle.ipa_tank_height = 0.85;
vehicle.n2o_tank_bottom = 0.5; % Bottom position of tanks
vehicle.ipa_tank_bottom = 0.5;
vehicle.cg_dry = 2; % Dry CG of vehicle
vehicle.mass_d = 50; % Dry mass (kg)
vehicle.Iyy_dry = 2500; % Dry moment of inertia (kg*m^2)
vehicle.tvc_pos = 0.1; % Vertical distance of TVC to ground
vehicle.initial_actuator_length = 0.17; % Neutral length of actuator
vehicle.actuator_gain = 4.5; % How much the TVC angle changes with the actuator length (rad/m)
vehicle.ure = 0; % Mass flow reference velocities
vehicle.wre = 0;

%% Inertia calculation: 
[Iyy, ~, ~, ~, ~, ~] = hopper_3dof_inertia_calculator(vehicle.n2o_mass_0, vehicle.ipa_mass_0, 0, vehicle);


%% Dynamics:
x = [0; 0; 0; 0; 0; 0; vehicle.n2o_mass_0; vehicle.ipa_mass_0; Iyy]; % Initial state
data = zeros(length(x)+8,N); % Store state history

% Simulation Loop
for i = 1:N
    u_ctrl = hopper_3dof_controller(x,i);
    [x, ctrl_data] = hopper_3dof_dynamics(x, u_ctrl, dt, vehicle);
    data(:,i) = [x ; ctrl_data];
end
t_vector = 0:dt:T;


hopper_3dof_plot_results(t_vector,data);


%% Trimming points:

u_val = zeros(2,N);
x_val = zeros(6,N);

for j = 1:N
uj= data( (16:17) ,  j ); % Thrust and deflection angle from "data" to lineraize about
xj = data( (1:6) ,  j ); % States from "data"

u_val(: , j) = uj;
x_val(: , j) = xj;
end 

 % loop over all the trimming points 
A_num = zeros(6,6,N);
B_num = zeros(6,2,N);
K     = zeros(2,6,N);
S     = zeros(6,6,N);
P     = zeros(6,1,N);

for i = 1:N
%% Get the Jacobian of the trimmed points: 
    [Ai, Bi] = HopperLinearization_3DOF(x_val(: , i), u_val(: , i), data(:,i));

    A_num(:,:,i) = Ai;
    B_num(:,:,i) = Bi;

    Q = eye(6);
    R = 0.1;
 %%  LQR controller 
    [K(:,:,i), S(:,:,i), P(:,:,i)] = lqr(Ai, Bi, Q, R);
end
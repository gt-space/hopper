
% Simulation Parameters
dt = 0.01; % Time step (s)
T = 10; % Time Span (s)
N = T/dt + 1; % Number of simulated states

% Vehicle Parameters
params.g = 9.81; % Gravity (m/s^2)

params.Tmax = 1556.88; % max thrust [N]
params.Isp = 260; % seconds
params.mass_w = 100; % Wet mass (kg)
params.mass_d = 50; % Dry mass (kg)

params.tvc_cg_dist_full = 1.2; % distance from CG to nozzle when full [m]
params.tvc_cg_dist_empty = 0.9; % distance from CG to nozzle when empty [m]
params.initial_actuator_length = 0.17; % Neutral length of actuator
params.actuator_gain = 4.5; % How much the TVC angle changes with the actuator length (rad/m)

params.ure = 0; % Mass flow reference velocities
params.wre = 0;

m = 100; % Mass (kg)
Iyy = 5000; % Moment of Inertia (kg*m^2)

x = [0; 0; 0; 0; 0; 0; m; Iyy]; % Initial state


data = zeros(length(x)+7,N); % Store state history

% Simulation Loop
for i = 1:N
    u_ctrl = hopper_3dof_controller(x);
    [x, ctrl_data] = hopper_3dof_dynamics(x, u_ctrl, dt, params);
    data(:,i) = [x ; ctrl_data];
end

t_vector = 0:dt:T;

hopper_3dof_plot_results(t_vector,data);




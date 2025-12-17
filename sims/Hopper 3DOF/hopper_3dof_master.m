
% Simulation Parameters
dt = 0.01; % Time step (s)
T = 20; % Time Span (s)
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

vehicle.tvc_pos = 0.1; % Vertical distance ofT TVC to ground
vehicle.initial_actuator_length = 0.17; % Neutral length of actuator
vehicle.actuator_gain = 4.5; % How much the TVC angle changes with the actuator length (rad/m)

vehicle.ure = 0; % Mass flow reference velocities
vehicle.wre = 0;

[Iyy, ~, ~, ~, ~, ~] = hopper_3dof_inertia_calculator(vehicle.n2o_mass_0, vehicle.ipa_mass_0, 0, vehicle);

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




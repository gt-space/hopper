function [x_next, ctrl_data] = hopper_3dof_dynamics(x, u_ctrl, dt, vehicle)

% Discrete-time 3DOF variable-mass body-axis plant
% Control inputs: throttle, actuator length

% --- Unpack state ---
u        = x(1); % X Velocity (m/s)
w        = x(2); % Z Velocity (m/s)
q        = x(3); % Pitch Angular Velocity (rad/s)
theta    = x(4); % Pitch Angle (rad)
xe       = x(5); % Inertial X Position (m)
ze       = x(6); % Inertial Z Position (m)
n2o_mass = x(7); % Mass (kg)
ipa_mass = x(8);
Iyy      = x(9); % Moment of Inertia (kg*m^2)

% --- Control inputs ---
throttle = max(0, min(1, u_ctrl.throttle)); % Throttle Percentage
actuator_length = vehicle.initial_actuator_length - u_ctrl.actuator_length; % Actuator Length (m)
phi = actuator_length * vehicle.actuator_gain; % Conversion to TVC Angle


% --- Thrust model ---
if n2o_mass >= 0 || ipa_mass >= 0
    T = throttle * vehicle.Tmax;
else
    T = 0; % Set thrust to zero if propellant is fully expended
end

% --- Mass and Inertia Calculations
mdot = -T / (vehicle.Isp * vehicle.g);
[Iyy, Iyy_dot, cg, mdot_n2o, mdot_ipa, m] = hopper_3dof_inertia_calculator(ipa_mass, n2o_mass, mdot, vehicle);

% --- TVC Calculations ---
lever = cg - vehicle.tvc_pos;
Tx =  T * sin(phi); % phi = 0 is TVC no gimbaled, + is to right
Tz = -T * cos(phi);
My =  lever * T * sin(phi);

% --- Continuous EOM ---
u_dot = Tx/m - q*w - (mdot/m)*(u - vehicle.ure) - vehicle.g*sin(theta);
w_dot = Tz/m + q*u - (mdot/m)*(w - vehicle.wre) + vehicle.g*cos(theta);
q_dot = (My - Iyy_dot*q)/Iyy;
theta_dot = q;

xe_dot =  u*cos(theta) - w*sin(theta);
ze_dot =  u*sin(theta) + w*cos(theta);

% --- Integrate (Euler) ---
x_next = x + dt * [
    u_dot;
    w_dot;
    q_dot;
    theta_dot;
    xe_dot;
    ze_dot;
    mdot_n2o;
    mdot_ipa;
    Iyy_dot
];


% Set Propellant masses to zero if fully expended
if (x(7) + dt * mdot_n2o) <= 0
    x_next(7) = 0;
end

if (x(8) + dt * mdot_ipa) <= 0
    x_next(8) = 0;
end

ctrl_data = [u_dot ; w_dot; q_dot; mdot; m; cg; T; phi];

end

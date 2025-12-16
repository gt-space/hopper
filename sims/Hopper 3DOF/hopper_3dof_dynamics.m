function [x_next, ctrl_data] = hopper_3dof_dynamics(x, u_ctrl, dt, params)

% Discrete-time 3DOF variable-mass body-axis plant
% Control inputs: throttle, actuator length

% --- Unpack state ---
u     = x(1); % X Velocity (m/s)
w     = x(2); % Z Velocity (m/s)
q     = x(3); % Pitch Angular Velocity (rad/s)
theta = x(4); % Pitch Angle (rad)
xe    = x(5); % Inertial X Position (m)
ze    = x(6); % Inertial Z Position (m)
m     = x(7); % Mass (kg)
Iyy   = x(8); % Moment of Inertia (kg*m^2)

% --- Control inputs ---
throttle = max(0, min(1, u_ctrl.throttle)); % Throttle Percentage
actuator_length = params.initial_actuator_length - u_ctrl.actuator_length; % Actuator Length (m)
phi = actuator_length * params.actuator_gain; % Conversion to TVC Angle


% --- TVC Lever Arm Calculation ---
lever = params.tvc_cg_dist_full + (params.tvc_cg_dist_empty - params.tvc_cg_dist_full) * ( 1 - m / (params.mass_w));

% --- Thrust model ---
if m > params.mass_d % Check if propellant is not fully expended
    T = throttle * params.Tmax;
    Tx =  T * sin(phi); % phi = 0 is TVC no gimbaled, + is to right
    Tz = -T * cos(phi);
    My =  lever * T * sin(phi);
else
    Tx = 0;
    Tz = 0;
    My = 0;
end

% --- Mass & inertia rates ---
mdot = -T / (params.Isp * params.g);
Iyy_dot = mdot * (lever^2);


% --- Continuous EOM ---
u_dot = Tx/m - q*w - (mdot/m)*(u - params.ure) - params.g*sin(theta);
w_dot = Tz/m + q*u - (mdot/m)*(w - params.wre) + params.g*cos(theta);
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
    mdot;
    Iyy_dot
];


if (x(7) + dt * mdot) <= params.mass_d
    x_next(7) = params.mass_d;
end

ctrl_data = [u_dot ; w_dot; q_dot; mdot; lever; T; phi];

end

function hopper_3dof_plot_results(t, X)

% Plot state history

u     = X(1,:);
w     = X(2,:);
q     = X(3,:);
theta = X(4,:);
xe    = X(5,:);
ze    = X(6,:);
m     = X(7,:);
Iyy   = X(8,:);
u_dot = X(9,:);
w_dot = X(10,:);
q_dot = X(11,:);
m_dot = X(12,:);
lever =  X(13,:);
thrust = X(14,:);
phi = X(15,:);


figure('Name','Hopper 3DOF Rotational & Translation States','Color','k');

%--Positions--
% --- Positions ---
subplot(3,3,1)
plot(t, xe, 'LineWidth',1.5)
grid on
ylabel('X Position (m)')
title('Downrange Position')

subplot(3,3,2)
plot(t, ze, 'LineWidth',1.5)
grid on
ylabel('Z Position (m)')
title('Altitude Position')

% --- Velocities ---
subplot(3,3,4)
plot(t, u, 'LineWidth',1.5)
grid on
ylabel('u (m/s)')
title('Body X Velocity')

subplot(3,3,5)
plot(t, -w, 'LineWidth',1.5)
grid on
ylabel('w (m/s)')
title('Body Z Velocity')

% --- Accelerations ---
subplot(3,3,7)
plot(t, u_dot, 'LineWidth',1.5)
grid on
ylabel('u_{dot} (m/s^2)')
title('Body X Acceleration')

subplot(3,3,8)
plot(t, -w_dot, 'LineWidth',1.5)
grid on
ylabel('w_{dot} (m/s^2)')
title('Body Z Acceleration')

% --- Attitude ---
subplot(3,3,3)
plot(t, rad2deg(theta), 'LineWidth',1.5)
grid on
ylabel('\theta (deg)')
title('Pitch Angle')

subplot(3,3,6)
plot(t, rad2deg(q), 'LineWidth',1.5)
grid on
ylabel('q (deg/s)')
title('Pitch Rate')
subplot(3,3,9)
plot(t, q_dot, 'LineWidth',1.5)
grid on
ylabel('q_{dot} (deg/s^2)')
title('Pitch Acceleration')

figure('Name','Hopper 3DOF Trajectory','Color','k');

plot(xe, -ze, 'LineWidth',1.5)
grid on
xlabel('Downrange X (m)')
ylabel('Altitude Z (m)')
title('Trajectory')

figure('Name','Hopper 3DOF Mass States','Color','k');
% --- Mass ---
subplot(2,2,1)
plot(t, m, 'LineWidth',1.5)
grid on
ylabel('Mass (kg)')
xlabel('Time (s)')
title('Mass Depletion')

% --- Mass Rate ---
subplot(2,2,2)
plot(t, m_dot, 'LineWidth',1.5)
grid on
ylabel('m_{dot} (kg/s)')
xlabel('Time (s)')
title('Mass Depletion Rate')

% --- Lever Arm ---
subplot(2,2,3)
plot(t, lever, 'LineWidth',1.5)
grid on
ylabel('Lever Arm (m)')
xlabel('Time (s)')
title('Lever Arm Profile')

% --- Moment of Inertia ---
subplot(2,2,4)
plot(t, Iyy, 'LineWidth',1.5)
grid on
ylabel('I_{yy} (kg*m^2)')
xlabel('Time (s)')
title('Moment of Inertia Profile')

figure('Name','Hopper 3DOF Control Inputs','Color','k');

% --- Thrust ---
subplot(1,2,1)
plot(t, thrust, 'LineWidth',1.5)
grid on
ylabel('Thrust (N)')
xlabel('Time (s)')
title('Thrust Profile')

% --- TVC Gimbal Angle ---
subplot(1,2,2)
plot(t, rad2deg(phi), 'LineWidth',1.5)
grid on
ylabel('\phi (deg)')
title('TVC Pitch Angle')

end
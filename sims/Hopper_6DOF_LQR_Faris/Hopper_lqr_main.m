%% Initial conditions:
% x(:,1:3)     = [X Y Z]              (position NED)
% x(:,4:6)     = [u v w]                (body velocity)
% x(:,7:9)     = [P Q R]                (body rates)
% x(:,10:13) = [q0 q1 q2 q3]      (quaternion)

%% some ideas:
%1- We can sweep over different values of Q and R and track the error to get
% lower value 
% 2- Add the mass change from the looking table, add rate limit to thrust
% and limit thrust. 
% 3- Test different trajectory 
% 4- Add ground to the sim
% 5- Linearize around cg step instead of time step
% still need the tvc req input to track 

%% ==== INITIAL CONDITIONS ====
clc, close all, clear all

x0 = zeros(13,1);
x0(12) = sqrt(2)/2;
x0(10) = sqrt(2)/2;

%% ==== TIME ====
tf = 25;                         % Final time
dt = 0.01;                      % Time step
N = 200;                       % Number of points 
t = linspace(0 , tf , N); % Time vector 
nx = 13;                        % Number of states 
nu = 4;                           % Number of inputs 
%% ==== CONSTANTS ====
m0  = 114;                    % initial mass [kg]
mdot = 0.5;                   % mass flow rate [kg/s]
mf = 50;
D_m = 0.25;
r_m = D_m/2;
L_m = 4;
g   = ones(N, 1) * 9.81;

%% ==== TIME-VARYING MASS AND INERTIA====
m =  max(m0 - mdot*t, mf)'; 

%% ==== TRAJECTORY AND CONTROL INPUT====
% Desired Trajectory:

% Preallocate arrays
Zref = zeros(length(t),1);
Wref = zeros(length(t),1);
Xref = zeros(length(t),1);
Uref = zeros(length(t),1);

for k = 1:length(t)
    
    if t(k) <= 10
        % Ascend (upward)
        az = -5;                         % acceleration in z
        Zref (k) = -5 * t(k);     % from 0 to -50 m
        Wref (k) = -5;

    elseif t(k) <= 15
        % Hover
        Zref (k) = -50;
        Wref (k) = 0;

    elseif t(k) <= 25
        % Descend
        Zref (k) = -50 + 5*(t(k)-15);
        Wref (k) = 5;

    else
        Zref(k) = 0;
        Wref(k) = 0;
    end

end

    % for i = (90: length(t))
    % ax = 1; % acceleration in x
    % Xref (i) =   0.5 * ax * t(i)^2;
    % Uref (i) = ax*t(i) ;
    % end 
        
x_ref = zeros(length(t),13);
x_ref(:, 1) = Xref;
x_ref(:, 5) = Uref;
x_ref(:,3) = Zref;
x_ref(:,4) = Wref;
x_ref(:,12) = sqrt(2)/2;
x_ref(:,10) = sqrt(2)/2;


figure;
plot(t, x_ref(:,3), 'r--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Zref [m]');
title('Reference Z Trajectory');
grid on;

figure;
plot(x_ref(:,1), x_ref(:,3), 'r--', 'LineWidth', 1.5);
xlabel('Xref [m]');
ylabel('Zref [m]');
title('Reference Trajectory');
grid on;

figure;
plot(t, x_ref(:,4), 'b--', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Wref [m/s]');
title('Reference Velocity in Z Direction');
grid on;

% Desired control input:
Zddot_ref = zeros(N,1);
Xddot_ref = zeros(N,1);
Tref = zeros(N,1);
delta_p_ref = zeros(N,1);

% for i = 2:length(t)-1
% 
%     Zddot_ref(i) = (Wref(i+1) - Wref(i-1)) / (t(i+1) - t(i-1));
%     Xddot_ref(i) = (Uref(i+1) - Uref(i-1)) / (t(i+1) - t(i-1));
% 
%     Fz = m(i) * (g(i) - Zddot_ref(i));
%     Fx = m(i) * Xddot_ref(i);
% 
%     Tref(i) = sqrt(Fz^2 + Fx^2);
% 
%     delta_p_ref(i) = atan2(Fx, Fz);
% 
% end


for i = 2:length(t)-1
    
    Zddot_ref(i) = (Wref(i+1) - Wref(i-1)) / (t(i+1) - t(i-1));
    T_unsat = m(i).*g(i) - m(i).*Zddot_ref(i);
    Tref(i) = min( max(T_unsat, 890), 2400 );

end

% Preallocate arrays
delta_y_ref = zeros(N,1);
F_rcs_ref = zeros(N,1);

u_ref = [Tref, delta_p_ref, delta_y_ref, F_rcs_ref ];

figure;
plot(t, Tref, 'k-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Reference Thrust [N]');
title('Reference Thrust over Time');
grid on;

figure;
plot(t, rad2deg(delta_p_ref), 'k-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Reference delta p [deg]');
title('Reference TVC deflection for pitch attitude over Time');
grid on;


%% ========= LINEARIZATION (JACOBIANS)  =========
% Preallocate arrays
u_val = zeros(nu,N);
x_val = zeros(nx,N);
m_value = zeros(1,N);
A_num = zeros(nx,nx,N);
B_num = zeros(nx,nu,N);

for j = 1:N
    u_val = u_ref(j ,:)';  % Take the desired input value at current time step
    x_val = x_ref(j ,:)';  % Take the desired state values at current time step
    m_value = m(j , :);   % Mass value at current time step
    g_value =g(j , :);

    [Ai, Bi] = HopperLinearization_6DOF(x_val, u_val, m_value , g_value);
    A_num(: , : , j) = Ai;
    B_num(: , : , j) = Bi;
end


%% ========= RICCATI ========= 
Q_riccati = diag([10, 10, 100000, 50, 50, 50, 100, 100,  100, 100, 100, 100, 100]); 
% Q Penalizes state deviations, higher Q values lead to faster, more aggressive control
R_riccati = diag([0.009, 50, 50, 50]);
% R penalizes control effort, higher R values result in more conservative
% conrtol 

% Backward time vector:
time_vec = tf : -dt : 0;  

% Preallocate arrays
S = zeros(nx, nx, N);
V = zeros(nx, N);
W = zeros(1, N);

% Terminal conditions (at t_f)
xd = x_ref(end,:)'; % Desired terminal states 
Sf = eye(13);
S(:,:,1) = Sf;  % First element is t_f
V(:,1) = -Sf * xd;
W(1) = 0.5 * xd' * Sf * xd;

% =========  Backward integration loop =========
for i = 1:N-1
    % Current backward time
    t_current = t(i);
    
    % Get A, B at this time
    Ai = A_num(:,:,i);
    Bi = B_num(:,:,i);
    
    % Get x_ref at this time
    xi_ref = x_ref(i, :)';
    
    % Get current S, V, W
    Si = S(:,:,i);
    Vi = V(:,i);
    Wi = W(i);
   
    % ----- Integrate S backward -----
    S_vec = Si(:);
    % Compute S derivative
    Sdot_vec = Riccati(t_current, S_vec, Ai, Bi, Q_riccati, R_riccati);
    
    S_next = reshape(S_vec - Sdot_vec * dt, [nx, nx]);
 
    S(:,:,i+1) = S_next;
    
    % ----- Integrate V backward -----
    Vdot = Riccati_V(t_current, Vi, Ai, Bi, Si, R_riccati, xi_ref);
    V(:,i+1) = Vi - Vdot * dt;
    
    % ----- Integrate W backward -----
    Wdot = Riccati_W(t_current, Wi, Bi, Q_riccati, R_riccati, xi_ref, Vi);
    W(i+1) = Wi - Wdot * dt;
end 

% Reverse arrays to get forward time
S = flip(S, 3);
V = flip(V, 2);
W = flip(W, 2);

%% ========= CONTROL CALCULATION =========
% Preallocate arrays
u_lqr = zeros(nu, N);
K = zeros (nu, nx , N);
 THRUST = zeros(N,1);
% Initialize states
x_sim = zeros(nx, N);
x_sim(:,1) = x0;  % Initial actual state:

for i = 1:N-1
    
    % Get current state from simulation:
    x_current_actual = x_sim(:, i);
    
    % Get current gains (S, V, B) at time t(i):
    Si = S(:,:,i);
    Vi = V(:,i);
    Bi = B_num(:,:,i);
    
    
    x_ref_i = x_ref(i,:)';
    u_ref_i = u_ref(i,:)';

    % Tracking error
    e = x_current_actual - x_ref_i;

    % Calculate gains at this time step:
    K(:, : , i) = inv(R_riccati) * Bi' * Si;

    %  tracking LQR control
    u_current = u_ref_i - K(:, : , i) * e - inv(R_riccati) * Bi' * Vi;

    %  Apply to 6DOF 
    T_current = min( max(u_current(1), 890), 2400 ); % 890 N to 2400 limit
    delta_p_current = u_current(2);
    delta_y_current = u_current(3);
    F_rcs_current = u_current(4);
    THRUST(i) = T_current;
    odefun = @(tt, xx) hopper6dof2(tt, xx, T_current, delta_p_current, ...
                                   delta_y_current, F_rcs_current, m(i), g(i));
    
    [~, x_temp] = ode45(odefun, [t(i) t(i+1)], x_current_actual);
    x_sim(:,i+1) = x_temp(end,:)';
    
    u_lqr(:,i) = u_current;
end

u_lqr(:,N) = u_lqr(:,N-1);


%% ============ PLOTS ============
figure
plot(t, -x_sim(3,:), 'b-', 'LineWidth', 2); hold on;
plot(t, -x_ref(:,3), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Altitude Tracking');
legend('Actual', 'Desired', 'Location', 'best');
grid on;


figure('Position', [100, 100, 1400, 900]);

subplot(3,3,1);
plot(t, -x_sim(3,:), 'b-', 'LineWidth', 2); hold on;
plot(t, -x_ref(:,3), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Altitude Tracking');
legend('Actual', 'Desired', 'Location', 'best');
grid on;

subplot(3,3,2);
plot(t, x_sim(4,:), 'b-', 'LineWidth', 2); hold on;
plot(t, x_ref(:,4), 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Vertical Velocity (m/s)');
title('Vertical Velocity (u in body frame)');
legend('Actual', 'Desired');
grid on;

subplot(3,3,3);
plot(t(1:end-1)', THRUST(1:end-1,1), 'k-', 'LineWidth', 2);
hold on;
plot([t(1) t(end)], [m(1)*9.81 m(1)*9.81], 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Thrust (N)');
title('Thrust Command');
legend('Command', 'Hover Thrust');
grid on;

euler_angles = quat2eul(x_sim(10:13,:)');
subplot(3,3,4);
plot(t, rad2deg(euler_angles(:,1)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, rad2deg(euler_angles(:,2)), 'r-', 'LineWidth', 1.5);
plot(t, rad2deg(euler_angles(:,3)), 'g-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (deg)');
title('Euler Angles');
legend('Roll (\phi)', 'Pitch (\theta)', 'Yaw (\psi)');
grid on;

subplot(3,3,5);
plot(t, (x_sim(7,:)), 'b-', 'LineWidth', 1.5); hold on;
plot(t, (x_sim(8,:)), 'r-', 'LineWidth', 1.5);
plot(t, (x_sim(9,:)), 'g-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Rate (deg/s)');
title('Body Rates');
legend('P', 'Q', 'R');
grid on;


subplot(3,3,6);
plot(t(1:end-1), (u_lqr(2,1:end-1)), 'b-', 'LineWidth', 2); hold on;
plot(t(1:end-1), (u_lqr(3,1:end-1)), 'r-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Command (deg)');
title('Pitch/Yaw Commands');
legend('\delta_p', '\delta_y');
grid on;

subplot(3,3,7);
alt_error = -x_sim(3,:) - (-x_ref(:,3)');
plot(t, alt_error, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Error (m)');
title('Altitude Tracking Error');
grid on;

subplot(3,3,8);
plot(t, x_sim(5,:), 'b-', 'LineWidth', 1.5); hold on;
plot(t, x_sim(6,:), 'r', 'LineWidth', 1.5);
hold on;
xlabel('Time (s)'); ylabel('Velocities');
title('v and w velocities');
grid on;

subplot(3,3,9);
quat_norm = sqrt(sum(x_sim(10:13,:).^2, 1));
plot(t, quat_norm, 'b-', 'LineWidth', 2);
hold on;
plot([t(1) t(end)], [1 1], 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('‖q‖');
title('Quaternion Norm');

grid on;


figure('Position', [100, 100, 800, 600]);
plot3(x_sim(1,:), x_sim(2,:), -x_sim(3,:), 'b-', 'LineWidth', 2);
hold on;
plot3(x_ref(:,1), x_ref(:,2), -x_ref(:,3), 'r--', 'LineWidth', 1.5);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Altitude (m)');
title('3D Trajectory');
legend('Actual', 'Desired');
grid on; view(45, 30);

%% Initial conditions:
% x(:,1:3)     = [X Y Z]              (position NED)
% x(:,4:6)     = [u v w]                (body velocity)
% x(:,7:9)     = [P Q R]                (body rates)
% x(:,10:13) = [q0 q1 q2 q3]      (quaternion)
clc, close all, clear all

x0 = zeros(13,1);
x0(12) = 0.7071;
x0(10) = 0.7071;

%% Define the trajectory:
Time = 25;
dt = 0.5;
N = Time/dt;
t = linspace(0, Time, N);
%% ==== CONSTANTS ====
m0  = 100;          % initial mass [kg]
mdot = 0.5;         % mass burn rate [kg/s]

D_m = 0.25;
r_m = D_m/2;
L_m = 4;
g   = ones(N, 1) * 9.81;

%% ==== TIME-VARYING MASS ====
m = max(m0 - mdot*t, 1)'; 

% Time vector from simulation
Zref = zeros(length(t),1);
Wref = zeros(length(t),1);

for k = 1:length(t)
    if t(k) <= 10
        % Ascend (upward)
        Zref(k) = -5 * t(k);     % from 0 to -50 m
        Wref(k) = -5;

    elseif t(k) <= 15
        % Hover
        Zref(k) = -50;
        Wref(k) = 0;

    elseif t(k) <= 25
        % Descend
        Zref(k) = -50 + 5*(t(k)-15);
        Wref(k) = 5;

    else
        Zref(k) = 0;
        Wref(k) = 0;
    end
end

x_ref = zeros(length(t),13);

% Position
x_ref(:,3) = Zref;

% Velocity (body frame)
x_ref(:,6) = Wref;

% Quaternion (level hover)
x_ref(:,12) = 0.7071;
x_ref(:,10) = 0.7071;


T = 1000;
delta_p = 0;
delta_y = 0;
F_rcs = 0;

%% ODE to solve the dynamics
[t, X] = ode45(@(t,x) hopper6dof(t, x, T, delta_p, delta_y, F_rcs, m0, r_m, L_m, g, mdot, m), t, x0);

%% Compute errors between desired trajectory and nominal flight:
x_err = X - x_ref;

%% Calculate Thrust profile based on trajectory 
Zddot_ref = zeros(size(t));

for k = 2:length(t)-1
    Zddot_ref(k) = (Wref(k+1) - Wref(k-1)) / (t(k+1) - t(k-1));
end

Tref = m .* g - m .* Zddot_ref; % Thrust required to achieve the trajectory (based on the trajectory)
delta_p_ref = zeros(N,1);
delta_y_ref = zeros(N,1);
F_rcs_ref = zeros(N,1);

figure;
plot(t, Tref, 'LineWidth', 1.5); hold on;
grid on;
xlabel('Time [s]')
ylabel('T [N]')
title('Desired thrust profile (N)')

u = [Tref, delta_p_ref, delta_y_ref, F_rcs_ref ]; % All desired inputs 
%% ========= JACOBIANS =========
% Trimming points:
u_val = zeros(4,N);
x_val = zeros(13,N);
m_value = zeros(1,N);

for j = 1:N
uj= u( j ,  : ); % Inputs
xj = x_ref( j , : ); % States 
mj = m(j, :); % mass

u_val(: , j) = uj;
x_val(: , j) = xj;
m_value(: , j) = mj;
end 


 % loop over all trimming points 
A_num       = zeros(13,13,N);
B_num        = zeros(13,4,N);

%% Get the Jacobian of the trimmed points: 

for i = 1:N
    [Ai, Bi] = HopperLinearization_6DOF(x_val(: , i), u_val(: , i) , m_value(: , i));

    A_num(:,:,i) = Ai;
    B_num(:,:,i) = Bi;
end

%%  ========= TIME-VARYING LQR (TV-LQR) =========
Q_riccati =  diag([0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);   % A in Bryson book, state weights 
R_riccati = diag([2, 1, 1, 1]);                                                                           % B in Bryson book, input weights

% Solving Riccati backward in time:
tf = t(length(t), :) : -dt : 0;

for s = N : -1 : 1
    % Get current A and B matrices:
    A = A_num(: , : , s);
    B = B_num(: , : , s);
 
    S_tf = eye(13); % Initial (Last) condition for the Riccati equation
    [tS, Sdot] = ode45(@(t,S) Riccati(tf, S, A, B, Q_riccati, R_riccati), tf, S_tf); % Backward integration for Riccati equation
end 

    S = reshape((Sdot).', [13, 13, 51]); % Reshaping into 13,13 , 51 (ode45 take vector only),

    % Solving for R:
 for r =  N : -1 : 1
     % Get current A, B and S matrices:
    Ai = A_num(:,:,r);
    Bi = B_num(:,:,r);
    S_current = S(:,:,r);

    R_tf = eye(13); % Initial (Last) condition 
    [tR, Rdot] = ode45(@(t,R) R_matrix(tf, S_current, R, Ai, Bi, Q_riccati, R_riccati ), tf, R_tf); % Backward integration to get R

end 
    R = reshape((Rdot).', [13, 13, 51]); % Reshaping into 13,13 , 51 (ode45 take vector only),

    % Solving for Q:
 for q =  N : -1 : 1
     % Get current A, B and S matrices:
    Ai = A_num(:,:,q);
    Bi = B_num(:,:,q);
    R_current = R(:,:,q);

    Q_tf = zeros(13); % Initial (Last) condition 
    [tQ, Qdot] = ode45(@(t,Q) Q_matrix(tf, R_current, Ai, Bi, Q_riccati, R_riccati ), tf, Q_tf); % Backward integration to get Q

end 
    Q = reshape((Qdot).', [13, 13, 51]); % Reshaping into 13,13 , 51 (ode45 take vector only),

    K= zeros(4,13,N);
 for i =  N : -1 : 1
     % Get current A, B and S matrices:
    Ai = A_num(:,:,i);
    Bi = B_num(:,:,i);
    S_current = S(:,:,i);
    R_current = R(:,:,i);
    Q_current = Q(:,:,i);

    Ki = pinv(R_riccati) * (Bi' * (S_current - R_current * (pinv(Q_current) * R_current')));
    K(:,:,i) = Ki;
 end 
 
K_flipped = K(:, :, end:-1:1);  % flip pages
states = reshape(X.', [13, 1, 50]);
X_error = reshape(x_err.', [13, 1, 50]);

u_closed = zeros(4, 1, 50);     % preallocate
for i = 1:50
    u_closed(:, :, i) = -K_flipped(:, :, i) * X_error(:, :, i);
end


%% Solve the dynamics for the closed loop
% Closed-loop simulation using TV-LQR 
[t, X_cl] = ode45(@(ti, x) ...
    hopper6dof(ti, x, ...
        u_closed(1,1,max(1,min(N,round(ti/dt)+1))), ...
        u_closed(2,1,max(1,min(N,round(ti/dt)+1))), ...
        u_closed(3,1,max(1,min(N,round(ti/dt)+1))), ...
        u_closed(4,1,max(1,min(N,round(ti/dt)+1))), ...
        m0, r_m, L_m, g(max(1,min(N,round(ti/dt)+1))), mdot, m(max(1,min(N,round(ti/dt)+1))) ...
    ), t, x0);


figure;
plot(t, x_ref(:,3), 'r--', 'LineWidth', 1.5); hold on;   % Reference Z
plot(t, X_cl(:,3), 'b', 'LineWidth', 1.5);               % Closed-loop Z
xlabel('Time [s]');
ylabel('Z Altitude [m]');
legend('Reference', 'Closed-loop');
grid on;
title('Z Altitude Tracking');
%%  ========= PLOTS =========
figure;
subplot(3,1,1)
plot(t, X(:,1)); grid on;
ylabel('X [m]')

subplot(3,1,2)
plot(t, X(:,2)); grid on;
ylabel('Y [m]')

subplot(3,1,3)
plot(t, X(:,3)); grid on;
ylabel('Z [m]')
xlabel('Time [s]')
sgtitle('Position (NED)')

figure;
subplot(3,1,1)
plot(t, X(:,4)); grid on;
ylabel('u [m/s]')

subplot(3,1,2)
plot(t, X(:,5)); grid on;
ylabel('v [m/s]')

subplot(3,1,3)
plot(t, X(:,6)); grid on;
ylabel('w [m/s]')
xlabel('Time [s]')
sgtitle('Body Velocities')

figure;
subplot(3,1,1)
plot(t, X(:,7)); grid on;
ylabel('P [rad/s]')

subplot(3,1,2)
plot(t, X(:,8)); grid on;
ylabel('Q [rad/s]')

subplot(3,1,3)
plot(t, X(:,9)); grid on;
ylabel('R [rad/s]')
xlabel('Time [s]')
sgtitle('Body Angular Rates')

figure;
plot(t, X(:,10), 'LineWidth', 1.5); hold on;
plot(t, X(:,11), 'LineWidth', 1.5);
plot(t, X(:,12), 'LineWidth', 1.5);
plot(t, X(:,13), 'LineWidth', 1.5);
grid on;

xlabel('Time [s]')
ylabel('Quaternion')
legend('q0','q1','q2','q3')
title('Quaternion States')

qnorm = sqrt(sum(X(:,10:13).^2, 2));

figure;
plot(t, qnorm, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]')
ylabel('||q||')
title('Quaternion Norm (should be 1)')

figure;
subplot(3,1,1)
plot(t, x_err(:,1)); grid on;
ylabel('e_X [m]')

subplot(3,1,2)
plot(t, x_err(:,2)); grid on;
ylabel('e_Y [m]')

subplot(3,1,3)
plot(t, x_err(:,3)); grid on;
ylabel('e_Z [m]')
xlabel('Time [s]')
sgtitle('Position Tracking Error')

figure;
plot(t, x_err(:,6), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]')
ylabel('e_w [m/s]')
title('Vertical Velocity Error')

figure;
plot(t, X(:,3), 'LineWidth', 1.5); hold on;
plot(t, Zref, '--', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]')
ylabel('Z [m]')
legend('Actual','Desired')
title('Altitude Tracking')

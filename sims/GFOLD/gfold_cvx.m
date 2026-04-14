%% GFOLD - Fuel Optimal Powered Descent Guidance
%  Convex optimization via CVX (SOCP formulation)
%
%  Based on: Acikmese & Ploen (2007), Blackmore et al. (2010)
%  "Lossless Convexification" approach - Problem 3 (min fuel)
%
% =========================================================
%  CVX SETUP (one-time, do this before running):
%
%  1. Download CVX from: http://cvxr.com/cvx/download/
%     -> Get "CVX for MATLAB" (free academic license)
%     -> Extract the zip anywhere, e.g. C:\cvx\
%
%  2. In MATLAB, navigate to that folder and run:
%        >> cd C:\cvx
%        >> cvx_setup
%     This installs CVX and its solvers (SDPT3, SeDuMi).
%
%  3. (Recommended) Get a free academic license at:
%        http://cvxr.com/cvx/academic/
%     Run >> cvx_setup('cvx_license.mat') with the file they email you.
%
%  4. Verify install: >> cvx_version
%     Should print solver info with no errors.
%
%  After that, this script runs with >> gfold_cvx
% =========================================================

clear; clc; close all;
tic;

%% ── Mission Parameters ───────────────────────────────────────────────────
% Source: mission_inputs.m + Hopper mass budget
% Phase: Ascent only (ground → 51m apex)

% Initial state
r0 = [0;  0;  0];        % position [m], z = altitude (start on ground)
v0 = [0;  0;  0];        % velocity [m/s] (start from rest)

% Terminal state
rf = [0;  0;  51];       % target position [m], z = 51m apex (mission_inputs.target_altitude)
vf = [0;  0;   0];       % velocity [m/s] (arrive at apex with zero velocity)

% Rate limits (mission_inputs.max_ascent_vel)
vz_max = 7;              % max vertical speed [m/s]

% Vehicle mass
m0      = 135.90;        % wet mass [kg]
m_prop  = 16.0;          % total propellant [kg] (ox=16 + fuel=16 per mission_inputs... 
                         % NOTE: ox_mass + fuel_mass = 32 kg in file but actual
                         % loaded prop = 16 kg per propulsion team)
m_dry   = m0 - m_prop;  % dry mass [kg] = 119.90 kg

% Propulsion (mission_inputs: nominal_thrust = 500 lbf, throttle_range = [0.4, 1.1])
T_nominal = 500 * 4.44822;          % N = 2224.1 N
Tmin      = 0.40 * T_nominal;       % N = 889.6 N  (40% throttle floor)
Tmax      = 1.10 * T_nominal;       % N = 2446.5 N (110% max)
Isp       = 225;                    % specific impulse [s] — placeholder, update when known
g0        = 9.80665;                % m/s^2 (mission_inputs.const.g0)
alpha     = 1 / (Isp * g0);        % mass flow coefficient [kg/(N·s)]

% Environment
g_vec = [0; 0; -9.80665];          % gravity vector [m/s^2]

% Constraints
theta_max = deg2rad(10);            % max thrust tilt from vertical — pure ascent, near 0
% Note: TVC max gimbal = 15 deg (mission_inputs.tvc.max_gimbal_angle)
% Using 10 deg here as conservative bound for trajectory planning

%% ── Discretization ───────────────────────────────────────────────────────

N  = 50;          % number of time steps
tf = 15;          % flight time [s] — 51m hop from rest, tuned for feasibility
dt = tf / N;      % time step [s]

% Time vector
t = linspace(0, tf, N+1);

%% ── Precompute Mass Bounds ───────────────────────────────────────────────
% z(k) = log(m(k)); mass bounded between dry mass and wet mass at each node

z0_low  = log(m0 - alpha * Tmax * t);          % max depletion rate
z0_up   = log(m0 - alpha * Tmin * t);          % min depletion rate
z0_low  = max(z0_low, log(m_dry));             % hard floor: never drop below dry mass

%% ── CVX Problem ──────────────────────────────────────────────────────────

fprintf('Setting up CVX problem (N=%d, tf=%.1f s)...\n', N, tf);

cvx_begin quiet
    % Decision variables
    variable r(3, N+1)      % position trajectory
    variable v(3, N+1)      % velocity trajectory
    variable u(3, N+1)      % thrust/mass = T/m (specific thrust) [m/s^2]
    variable z(1, N+1)      % log(mass) at each node
    variable sig(1, N+1)    % slack: ||u|| <= sig (lossless convexification)

    % ── Objective: maximize final mass = minimize fuel used ──
    maximize( z(N+1) )

    subject to

        % ── Initial conditions ──
        r(:,1) == r0;
        v(:,1) == v0;
        z(1)   == log(m0);

        % ── Terminal conditions ──
        r(:,N+1) == rf;
        v(:,N+1) == vf;

        % ── Dynamics (trapezoidal integration) ──
        for k = 1:N
            % Position
            r(:,k+1) == r(:,k) + dt/2 * (v(:,k) + v(:,k+1));
            % Velocity: v_dot = u + g
            v(:,k+1) == v(:,k) + dt/2 * ( (u(:,k) + g_vec) + (u(:,k+1) + g_vec) );
            % Log-mass: z_dot = -alpha * sig  (sig = ||T||/m, z=log m)
            z(k+1)   == z(k) - dt/2 * alpha * (sig(k) + sig(k+1));
        end

        % ── Thrust magnitude bounds (lossless convexification) ──
        % ||u|| <= sig  (SOC constraint, convex)
        for k = 1:N+1
            norm(u(:,k), 2) <= sig(k);
        end

        % Lower bound: Tmin/m <= sig
        % Upper bound: sig <= Tmax/m
        % In log-mass: Tmin*exp(-z) <= sig <= Tmax*exp(-z)
        % Linearized around z0_low for convexity (first-order Taylor):
        for k = 1:N+1
            z0k = z0_low(k);
            % sig >= Tmin * exp(-z0k) * (1 - (z(k) - z0k))   [lower]
            sig(k) >= Tmin * exp(-z0k) * (1 - (z(k) - z0k));
            % sig <= Tmax * exp(-z0k) * (1 - (z(k) - z0k))   [upper]
            sig(k) <= Tmax * exp(-z0k) * (1 - (z(k) - z0k));
        end

        % ── Log-mass bounds ──
        for k = 1:N+1
            z(k) >= log(m_dry);   % hard floor: mass cannot drop below dry mass
            z(k) <= z0_up(k);
        end

        % ── Pointing constraint: thrust within theta_max of vertical ──
        % u(3) >= ||u|| * cos(theta_max)   (z-component must dominate)
        for k = 1:N+1
            u(3,k) >= sig(k) * cos(theta_max);
        end

        % ── Vertical rate limit: |vz| <= vz_max ──
        for k = 1:N+1
            v(3,k) <= vz_max;
            v(3,k) >= -vz_max;
        end

        % ── No ground penetration ──
        for k = 1:N+1
            r(3,k) >= 0;
        end

        % ── sig must be positive ──
        sig >= 0;

cvx_end

%% ── Post-Process Results ─────────────────────────────────────────────────

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    fprintf('\nCVX Status: %s\n', cvx_status);

    % Recover mass and thrust from log-mass and specific thrust
    mass_traj   = exp(z);                        % [kg]
    thrust_traj = u .* mass_traj;                % T = u * m [N]
    thrust_mag  = sqrt(sum(thrust_traj.^2, 1));  % ||T|| [N]
    fuel_used   = m0 - mass_traj(end);

    printSummary(t, r, v, thrust_mag, mass_traj, m0, fuel_used, rf, vf);
    plotResults(t, r, v, thrust_traj, thrust_mag, mass_traj, m0, Tmax);

else
    fprintf('\n[!] CVX Status: %s\n', cvx_status);
    fprintf('Solver could not find a solution. Try:\n');
    fprintf('  1. Increase tf (current: %.1f s) — vehicle may need more time\n', tf);
    fprintf('  2. Decrease N — coarser grid sometimes helps feasibility\n');
    fprintf('  3. Relax theta_max or gamma_gs constraints\n');
    fprintf('  4. Check Tmin/Tmax ratio (Tmin=%.0f N, Tmax=%.0f N)\n', Tmin, Tmax);
end

totalTime = toc;
fprintf('Total computation time: %.3f s\n', totalTime);

%% ── Helper: Print Summary ────────────────────────────────────────────────
function printSummary(t, r, v, thrust_mag, mass_traj, m0, fuel_used, rf, vf)
    fprintf('\n============ GFOLD Solution Summary ============\n');
    fprintf('Flight time:       %.2f s\n', t(end));
    fprintf('Initial mass:      %.2f kg\n', m0);
    fprintf('Final mass:        %.2f kg\n', mass_traj(end));
    fprintf('Fuel used:         %.2f kg (%.1f%%)\n', fuel_used, 100*fuel_used/m0);
    fprintf('Max thrust:        %.1f N\n', max(thrust_mag));
    fprintf('Avg thrust:        %.1f N\n', mean(thrust_mag));
    fprintf('Initial position:  [%.2f, %.2f, %.2f] m\n', r(1,1), r(2,1), r(3,1));
    fprintf('Final position:    [%.2f, %.2f, %.2f] m\n', r(1,end), r(2,end), r(3,end));
    fprintf('Target position:   [%.2f, %.2f, %.2f] m\n', rf(1), rf(2), rf(3));
    fprintf('Final velocity:    [%.3f, %.3f, %.3f] m/s\n', v(1,end), v(2,end), v(3,end));
    fprintf('Target velocity:   [%.3f, %.3f, %.3f] m/s\n', vf(1), vf(2), vf(3));
    pos_err = norm(r(:,end) - rf);
    vel_err = norm(v(:,end) - vf);
    fprintf('Position error:    %.4f m\n', pos_err);
    fprintf('Velocity error:    %.4f m/s\n', vel_err);
    fprintf('================================================\n\n');
end

%% ── Helper: Plot Results ─────────────────────────────────────────────────
function plotResults(t, r, v, thrust_traj, thrust_mag, mass_traj, m0, Tmax)
    figure('Name','GFOLD Results','Position',[50 50 1400 900]);

    % ── 3D Trajectory ──
    subplot(2,3,1);
    plot3(r(1,:), r(2,:), r(3,:), 'b-', 'LineWidth', 2.5); hold on;
    plot3(r(1,1),   r(2,1),   r(3,1),   'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot3(r(1,end), r(2,end), r(3,end), 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % Draw thrust vectors (subsampled)
    skip = 3;
    scale = 0.004;
    for k = 1:skip:size(r,2)
        quiver3(r(1,k), r(2,k), r(3,k), ...
                thrust_traj(1,k)*scale, thrust_traj(2,k)*scale, thrust_traj(3,k)*scale, ...
                0, 'r', 'LineWidth', 1.2);
    end
    grid on; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('3D Trajectory + Thrust Vectors');
    legend('Trajectory','Start','Landing','Thrust','Location','best');
    view(45, 25);

    % ── Position vs Time ──
    subplot(2,3,2);
    plot(t, r(1,:), 'r-', t, r(2,:), 'g-', t, r(3,:), 'b-', 'LineWidth', 1.8);
    grid on; xlabel('Time [s]'); ylabel('Position [m]');
    legend('X','Y','Z (alt)','Location','best');
    title('Position vs Time');

    % ── Velocity vs Time ──
    subplot(2,3,3);
    plot(t, v(1,:), 'r-', t, v(2,:), 'g-', t, v(3,:), 'b-', 'LineWidth', 1.8);
    grid on; xlabel('Time [s]'); ylabel('Velocity [m/s]');
    legend('Vx','Vy','Vz','Location','best');
    title('Velocity vs Time');

    % ── Thrust Magnitude vs Time ──
    subplot(2,3,4);
    plot(t, thrust_mag, 'k-', 'LineWidth', 1.8); hold on;
    yline(Tmax, 'r--', sprintf('Tmax = %.0f N', Tmax), 'LineWidth', 1.2);
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust Magnitude vs Time');
    ylim([0, Tmax * 1.15]);

    % ── Mass vs Time ──
    subplot(2,3,5);
    plot(t, mass_traj, 'm-', 'LineWidth', 1.8); hold on;
    yline(m0, 'b--', sprintf('m0 = %.0f kg', m0), 'LineWidth', 1.2);
    grid on; xlabel('Time [s]'); ylabel('Mass [kg]');
    title('Mass vs Time');
    ylim([min(mass_traj)*0.95, m0*1.05]);

    % ── Thrust Components vs Time ──
    subplot(2,3,6);
    plot(t, thrust_traj(1,:), 'r-', ...
         t, thrust_traj(2,:), 'g-', ...
         t, thrust_traj(3,:), 'b-', 'LineWidth', 1.8);
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    legend('Tx','Ty','Tz','Location','best');
    title('Thrust Components vs Time');

    sgtitle('GFOLD: Fuel-Optimal Powered Descent', 'FontSize', 14, 'FontWeight', 'bold');
end
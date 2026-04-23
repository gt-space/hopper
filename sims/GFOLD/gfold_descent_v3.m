%% GFOLD - Fuel Optimal Powered Descent
%  Phase: DESCENT (51m apex -> 1m landing)
%
%  Formulation: Direct thrust vector [N], post-process mass
%  Key insight: short tf (8s) allows ballistic free-fall + braking
%  profile, which is both fuel-optimal and physically correct.
%  Long tf (20s) forces continuous thrusting to satisfy vz limit -> infeasible.
%
%  Source: mission_inputs.m + Hopper mass budget
%  Update ascent_fuel_used if ascent script (gfold_cvx.m) is re-run.
% =========================================================

clear; clc; close all;
tic;

%% ── Mission Parameters ───────────────────────────────────────────────────

% Ascent handoff
ascent_fuel_used = 8.76;             % kg — update if ascent re-run
m_wet_initial    = 135.90;           % kg

% States
r0 = [0; 0; 51];   v0 = [0; 0; 0];  % apex: at rest
rf = [0; 0;  0.05];   vf = [0; 0; 0];  % landing: 1m AGL, at rest

% Vehicle mass
m_prop_total  = 16.0;
m_dry         = m_wet_initial - m_prop_total;    % 119.90 kg
m0            = m_wet_initial - ascent_fuel_used; % 127.14 kg at apex
m_prop_remain = m0 - m_dry;                      % 7.24 kg for descent

% Propulsion (mission_inputs)
T_nominal = 500 * 4.44822;           % N
Tmin      = 0.40 * T_nominal;        % N = 778.4 N — 35% throttle floor
                                     % Note: 40% (890N) hits prop budget limit for this descent
                                     % 35% gives optimizer margin; engine stays well above idle
Tmax      = 1.10 * T_nominal;        % N = 2446.5 N
Isp       = 225;                     % s
g0        = 9.80665;                 % m/s^2

% Environment & constraints
g_vec     = [0; 0; -g0];
vz_max    = 7;                       % m/s (mission_inputs.max_descent_vel)
theta_max = deg2rad(15);             % TVC gimbal limit (mission_inputs.tvc.max_gimbal_angle)

%% ── Discretization ───────────────────────────────────────────────────────
% Short tf forces ballistic free-fall + braking burn profile
% This is the fuel-optimal shape and stays within prop budget

N  = 60;
tf = 12;   % s — min feasible tf with 35% Tmin is ~8.76s; 12s gives optimizer margin
dt = tf / N;
t  = linspace(0, tf, N+1);

fprintf('Propellant available for descent: %.2f kg\n', m_prop_remain);
fprintf('Tmin=%.0fN (%.0f%% throttle), Tmax=%.0fN\n', Tmin, 100*Tmin/T_nominal, Tmax);
fprintf('Setting up CVX problem (N=%d, tf=%.1f s)...\n', N, tf);

%% ── CVX Problem ──────────────────────────────────────────────────────────
% Decision variables: thrust vector T [N], magnitude slack Tm [N]
% Mass tracked post-solve from Tm profile (avoids SOC+equality chain issue)
% Kinematics linearized at m0 (conservative, keeps problem convex)

cvx_begin quiet
    variable r(3, N+1)      % position [m]
    variable v(3, N+1)      % velocity [m/s]
    variable T(3, N+1)      % thrust vector [N]
    variable Tm(1, N+1)     % thrust magnitude [N] (SOC slack)

    % Minimize total thrust impulse = minimize fuel
    minimize( sum(Tm) * dt )

    subject to

        % ── Boundary conditions ──
        r(:,1) == r0;    v(:,1) == v0;
        r(:,N+1) == rf;  v(:,N+1) == vf;

        % ── Dynamics (trapezoidal, linearized at m0) ──
        for k = 1:N
            r(:,k+1) == r(:,k) + dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1) == v(:,k) + dt/2*( ...
                (T(:,k)/m0 + g_vec) + (T(:,k+1)/m0 + g_vec));
        end

        % ── Thrust SOC ──
        for k = 1:N+1
            norm(T(:,k), 2) <= Tm(k);
        end

        % ── Thrust magnitude bounds ──
        Tm >= Tmin;   % 40% throttle floor — engine always on during descent
        Tm <= Tmax;

        % ── Pointing: thrust within theta_max of vertical ──
        for k = 1:N+1
            T(3,k) >= Tm(k) * cos(theta_max);
        end

        % ── Descent rate limit ──
        for k = 1:N+1
            v(3,k) >= -vz_max;
            v(3,k) <=  vz_max;
        end

        % ── Propellant budget (global linear constraint) ──
        sum(Tm) * dt <= m_prop_remain * Isp * g0;

        % ── No ground penetration ──
        r(3,:) >= 0;
        Tm     >= 0;

cvx_end

%% ── Post-Process ─────────────────────────────────────────────────────────

if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('\nCVX Status: %s\n', cvx_status);

    % Reconstruct mass from thrust profile
    m_traj = zeros(1, N+1);
    m_traj(1) = m0;
    for k = 1:N
        m_traj(k+1) = m_traj(k) - dt/(2*Isp*g0)*(Tm(k)+Tm(k+1));
    end

    thrust_mag = sqrt(sum(T.^2, 1));
    fuel_used  = m0 - m_traj(end);

    fprintf('\n======== GFOLD Descent Summary ========\n');
    fprintf('Flight time:            %.2f s\n', tf);
    fprintf('Initial mass:           %.2f kg\n', m0);
    fprintf('Final mass:             %.2f kg\n', m_traj(end));
    fprintf('Fuel used:              %.2f kg (%.1f%% of available)\n', ...
            fuel_used, 100*fuel_used/m_prop_remain);
    fprintf('Prop remaining (total): %.2f kg\n', m_prop_remain - fuel_used);
    fprintf('Max thrust:             %.1f N\n', max(thrust_mag));
    fprintf('Avg thrust:             %.1f N\n', mean(thrust_mag));
    fprintf('Max |vz|:               %.3f m/s\n', max(abs(v(3,:))));
    fprintf('Final position:         [%.3f, %.3f, %.3f] m\n', r(1,end),r(2,end),r(3,end));
    fprintf('Final velocity:         [%.3f, %.3f, %.3f] m/s\n', v(1,end),v(2,end),v(3,end));
    fprintf('Position error:         %.4f m\n', norm(r(:,end)-rf));
    fprintf('Velocity error:         %.4f m/s\n', norm(v(:,end)-vf));
    fprintf('========================================\n\n');

    plotDescent(t, r, v, T, thrust_mag, m_traj, m0, Tmax, vz_max);

else
    fprintf('\n[!] CVX Status: %s\n', cvx_status);
    fprintf('Try: decrease tf (now %.1fs) — too long forces continuous thrusting\n', tf);
    fprintf('     increase N (now %d) for finer resolution\n', N);
end

totalTime = toc;
fprintf('Computation time: %.3f s\n', totalTime);

%% ── Plot ─────────────────────────────────────────────────────────────────
function plotDescent(t, r, v, T, thrust_mag, m_traj, m0, Tmax, vz_max)
    figure('Name','GFOLD Descent','Position',[50 50 1400 900]);

    % 3D trajectory
    subplot(2,3,1);
    plot3(r(1,:), r(2,:), r(3,:), 'b-', 'LineWidth', 2.5); hold on;
    plot3(r(1,1),   r(2,1),   r(3,1),   'go', 'MarkerSize',12, 'MarkerFaceColor','g');
    plot3(r(1,end), r(2,end), r(3,end), 'rs', 'MarkerSize',12, 'MarkerFaceColor','r');
    skip=4; scale=0.003;
    for k=1:skip:size(r,2)
        quiver3(r(1,k),r(2,k),r(3,k), ...
                T(1,k)*scale, T(2,k)*scale, T(3,k)*scale, ...
                0,'r','LineWidth',1.2);
    end
    grid on; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('3D Trajectory + Thrust Vectors');
    legend('Trajectory','Apex','Landing','Thrust','Location','best');
    view(45,25);

    % Position vs time
    subplot(2,3,2);
    plot(t, r(1,:),'r-', t, r(2,:),'g-', t, r(3,:),'b-', 'LineWidth',1.8);
    yline(1,'k--','Landing (1m)','LineWidth',1);
    grid on; xlabel('Time [s]'); ylabel('Position [m]');
    legend('X','Y','Z (alt)'); title('Position vs Time');

    % Velocity vs time
    subplot(2,3,3);
    plot(t, v(1,:),'r-', t, v(2,:),'g-', t, v(3,:),'b-', 'LineWidth',1.8); hold on;
    yline(-vz_max,'k--',sprintf('vz limit = -%.0f m/s',vz_max),'LineWidth',1);
    grid on; xlabel('Time [s]'); ylabel('Velocity [m/s]');
    legend('Vx','Vy','Vz'); title('Velocity vs Time');

    % Thrust magnitude
    subplot(2,3,4);
    plot(t, thrust_mag,'k-','LineWidth',1.8); hold on;
    yline(Tmax,'r--',sprintf('Tmax=%.0fN',Tmax),'LineWidth',1);
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust Magnitude vs Time');
    ylim([0, Tmax*1.15]);

    % Mass vs time
    subplot(2,3,5);
    plot(t, m_traj,'m-','LineWidth',1.8); hold on;
    yline(m0,'b--',sprintf('m0=%.0fkg',m0),'LineWidth',1);
    grid on; xlabel('Time [s]'); ylabel('Mass [kg]');
    title('Mass vs Time (post-processed)');
    ylim([min(m_traj)*0.99, m0*1.01]);

    % Thrust components
    subplot(2,3,6);
    plot(t, T(1,:),'r-', t, T(2,:),'g-', t, T(3,:),'b-','LineWidth',1.8);
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    legend('Tx','Ty','Tz'); title('Thrust Components vs Time');

    sgtitle('GFOLD: Fuel-Optimal Descent Phase','FontSize',14,'FontWeight','bold');
end
%% GFOLD - Fuel Optimal Powered Descent
%  Phase: DESCENT (51m apex -> 1m landing)
%  Formulation: Direct thrust vector, linear mass tracking
%  No log-mass substitution — avoids SOC+equality chain numerical issues
%
%  Source: mission_inputs.m + Hopper mass budget
%  Run gfold_cvx.m (ascent) first; update ascent_fuel_used if re-run.
% =========================================================

clear; clc; close all;
tic;

%% ── Mission Parameters ───────────────────────────────────────────────────

% Ascent handoff (update if ascent script is re-run)
ascent_fuel_used = 8.76;            % kg
m_wet_initial    = 135.90;          % kg

% States
r0 = [0; 0; 51];    v0 = [0; 0; 0];   % apex: at rest
rf = [0; 0;  1];    vf = [0; 0; 0];   % landing: 1m AGL, at rest

% Vehicle
m_prop_total  = 16.0;
m_dry         = m_wet_initial - m_prop_total;   % 119.90 kg
m0            = m_wet_initial - ascent_fuel_used; % 127.14 kg at apex
m_prop_remain = m0 - m_dry;                     % 7.24 kg for descent

% Propulsion
T_nominal = 500 * 4.44822;          % N
Tmin      = 0.0  * T_nominal;       % N — allow coast (engine off)
Tmax      = 1.10 * T_nominal;       % N = 2446.5 N
Isp       = 225;                    % s
g0        = 9.80665;                % m/s^2
mdot_max  = Tmax / (Isp * g0);     % max mass flow rate [kg/s]

% Environment & constraints
g_vec     = [0; 0; -g0];
vz_max    = 7;                      % m/s (mission_inputs.max_descent_vel)
theta_max = deg2rad(15);            % thrust tilt limit — TVC gimbal limit

%% ── Discretization ───────────────────────────────────────────────────────

N  = 60;
tf = 20;
dt = tf / N;
t  = linspace(0, tf, N+1);

fprintf('Propellant available: %.2f kg\n', m_prop_remain);
fprintf('Setting up CVX problem (N=%d, tf=%.1f s)...\n', N, tf);

%% ── CVX Problem ──────────────────────────────────────────────────────────
% Formulation: thrust vector T [N] as direct decision variable
% Kinematics:  a = T/m + g  (use fixed m0 for convexity — conservative)
% Mass:        m(k+1) = m(k) - dt/2*(Tm(k)+Tm(k+1))/(Isp*g0)
% SOC:         norm(T) <= Tm  (Tm is thrust magnitude slack)
% This avoids the log-mass + SOC interaction that breaks SDPT3

cvx_begin quiet
    variable r(3, N+1)      % position [m]
    variable v(3, N+1)      % velocity [m/s]
    variable T(3, N+1)      % thrust vector [N]
    variable Tm(1, N+1)     % thrust magnitude [N] (SOC slack)
    variable m(1, N+1)      % mass [kg]

    % Minimize fuel = minimize total thrust impulse
    minimize( sum(Tm) * dt )

    subject to

        % ── Boundary conditions ──
        r(:,1) == r0;   v(:,1) == v0;   m(1) == m0;
        r(:,N+1) == rf; v(:,N+1) == vf;

        % ── Dynamics (trapezoidal) ──
        % Use m0 for T/m to keep problem linear (outer relaxation)
        for k = 1:N
            r(:,k+1) == r(:,k) + dt/2*(v(:,k)+v(:,k+1));
            v(:,k+1) == v(:,k) + dt/2*( ...
                (T(:,k)/m0 + g_vec) + (T(:,k+1)/m0 + g_vec));
            m(k+1) == m(k) - dt/(2*Isp*g0)*(Tm(k)+Tm(k+1));
        end

        % ── Thrust SOC constraint ──
        for k = 1:N+1
            norm(T(:,k), 2) <= Tm(k);
        end

        % ── Thrust magnitude bounds ──
        Tm >= Tmin;
        Tm <= Tmax;

        % ── Pointing: thrust must be within theta_max of vertical ──
        for k = 1:N+1
            T(3,k) >= Tm(k) * cos(theta_max);
        end

        % ── Descent rate limit ──
        for k = 1:N+1
            v(3,k) >= -vz_max;
            v(3,k) <=  vz_max;
        end

        % ── Mass bounds ──
        m >= m_dry;
        m <= m0;

        % ── No ground penetration ──
        r(3,:) >= 0;

        % ── Thrust non-negative (magnitude) ──
        Tm >= 0;

cvx_end

%% ── Results ──────────────────────────────────────────────────────────────

if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    fprintf('\nCVX Status: %s\n', cvx_status);

    thrust_mag = sqrt(sum(T.^2, 1));
    fuel_used  = m0 - m(end);

    % Print summary
    fprintf('\n======== GFOLD Descent Summary ========\n');
    fprintf('Flight time:        %.2f s\n', tf);
    fprintf('Initial mass:       %.2f kg\n', m0);
    fprintf('Final mass:         %.2f kg\n', m(end));
    fprintf('Fuel used:          %.2f kg (%.1f%% of remaining)\n', ...
            fuel_used, 100*fuel_used/m_prop_remain);
    fprintf('Max thrust:         %.1f N\n', max(thrust_mag));
    fprintf('Avg thrust:         %.1f N\n', mean(thrust_mag));
    fprintf('Final position:     [%.3f, %.3f, %.3f] m\n', r(1,end),r(2,end),r(3,end));
    fprintf('Final velocity:     [%.3f, %.3f, %.3f] m/s\n', v(1,end),v(2,end),v(3,end));
    fprintf('Position error:     %.4f m\n', norm(r(:,end)-rf));
    fprintf('Velocity error:     %.4f m/s\n', norm(v(:,end)-vf));
    fprintf('Max |vz|:           %.3f m/s\n', max(abs(v(3,:))));
    fprintf('=======================================\n\n');

    plotDescent(t, r, v, T, thrust_mag, m, m0, Tmax);
else
    fprintf('\n[!] CVX Status: %s\n', cvx_status);
    fprintf('Try: increase tf (now %.1f s) or increase N (now %d)\n', tf, N);
end

totalTime = toc;
fprintf('Computation time: %.3f s\n', totalTime);

%% ── Plot ─────────────────────────────────────────────────────────────────
function plotDescent(t, r, v, T, thrust_mag, m, m0, Tmax)
    figure('Name','GFOLD Descent','Position',[50 50 1400 900]);

    subplot(2,3,1);
    plot(r(3,:), 'b-', 'LineWidth', 2); grid on;
    xlabel('Node'); ylabel('Altitude [m]');
    title('Altitude vs Node');
    yline(1, 'r--', 'Landing (1m)');

    subplot(2,3,2);
    plot(t, r(1,:),'r-', t, r(2,:),'g-', t, r(3,:),'b-', 'LineWidth',1.8);
    grid on; xlabel('Time [s]'); ylabel('Position [m]');
    legend('X','Y','Z (alt)'); title('Position vs Time');

    subplot(2,3,3);
    plot(t, v(1,:),'r-', t, v(2,:),'g-', t, v(3,:),'b-', 'LineWidth',1.8);
    hold on; yline(-7,'k--','vz limit'); yline(7,'k--');
    grid on; xlabel('Time [s]'); ylabel('Velocity [m/s]');
    legend('Vx','Vy','Vz'); title('Velocity vs Time');

    subplot(2,3,4);
    plot(t, thrust_mag, 'k-', 'LineWidth', 1.8); hold on;
    yline(Tmax,'r--',sprintf('Tmax=%.0fN',Tmax));
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    title('Thrust Magnitude vs Time');
    ylim([0, Tmax*1.15]);

    subplot(2,3,5);
    plot(t, m, 'm-', 'LineWidth', 1.8); hold on;
    yline(m0,'b--',sprintf('m0=%.0fkg',m0));
    grid on; xlabel('Time [s]'); ylabel('Mass [kg]');
    title('Mass vs Time');

    subplot(2,3,6);
    plot(t, T(1,:),'r-', t, T(2,:),'g-', t, T(3,:),'b-', 'LineWidth',1.8);
    grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
    legend('Tx','Ty','Tz'); title('Thrust Components vs Time');

    sgtitle('GFOLD: Fuel-Optimal Descent Phase','FontSize',14,'FontWeight','bold');
end
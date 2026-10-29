% State order:
% [X Y Z Vx Vy Vz P Q R p1 p2 p3]'
%
% Input order:
% [T delta_p delta_y Tor_rcs]'
%
% Integral error convention:
% xi_dot = position_command - position

Qdiags  = [0 0 0 10 10 100 10 10 10 1e6 1e6 1e6];
QIdiags = [100 100 500];

Q  = diag(Qdiags);
QI = diag([Qdiags QIdiags]);

Rdiags = [1 1e8 1e8 1e8];
R = diag(Rdiags);

%% Parameters

params = struct();
params.g = 9.81;
params.d = cg_init - engine_cg;
params.r_m = 0.125;
params.n_rcs = 1;

% MoI_init must be the ORIGINAL inertia in the OLD body axes.
% Do not apply this transformation twice.
S = [ 0  0  1;
      0  1  0;
     -1  0  0];

params.I = S * MoI_init * S.';

% Use this same transformed inertia in the Simulink 6DOF block.
I_sim = params.I;

m_profile = linspace(OUT.Vehicle.WetMass, ...
                     OUT.Vehicle.DryMass, length(t));
m = m_profile(1);

%% Hover equilibrium

p = zeros(3,1);

xi = [zeros(9,1); p];
ui = [m*params.g; 0; 0; 0];

[A,B] = HopperLinearization_lqi(xi,ui,params,m);

%% Integral augmentation: xi_dot = r - y

Cpos = [eye(3), zeros(3,9)];

Atilde = [ A,     zeros(12,3);
          -Cpos,  zeros(3,3)];

Btilde = [B; zeros(3,4)];

% Aliases for your linear-response code.
Aa = Atilde;
Ba = Btilde;

CrossTerm = zeros(15,4);

% RAW MATLAB gains:
% u = ubar - Klqi*[x-xbar; xi]
[KlqI,~,~] = lqr(Atilde,Btilde,QI,R,CrossTerm)


%% Physical Plant and LQI Controller for Disk Margin Analysis

% Extract physical-state and integrator gains
Kx = KlqI(:,1:12);    % 4 x 12
Ki = KlqI(:,13:15);   % 4 x 3

%% Plant P
% 4 actuator inputs -> 12 physical states

P = ss(A, B, eye(12), zeros(12,4));

%% Controller C
% Integral dynamics: xi_dot = -Cpos*x
%
% Controller output v = Kx*x + Ki*xi
% Negative feedback applies u = -v.

Ac = zeros(3,3);
Bc = -Cpos;
Cc = Ki;
Dc = Kx;

C = ss(Ac, Bc, Cc, Dc);

%% Verify the closed-loop dynamics

Acl_expected = Atilde - Btilde*KlqI;

% MATLAB's feedback uses negative feedback by default
CL = feedback(P,C);

poles_expected = sort(eig(Acl_expected));
poles_actual = sort(eig(CL.A));

% Compare characteristic polynomials via eigenvalue sets
poleError = max(abs(poles_expected - poles_actual));

fprintf('Closed-loop pole error: %.3e\n', poleError);


%% Disk Margin Analysis

sigma = 0; %skew

% Plant output margins (12 channels)
[DMo, MMo] = diskmargin(P*C, sigma);

% Plant input margins (4 channels)
[DMi, MMi] = diskmargin(C*P, sigma);

% Simultaneous input and output margins (16 channels)
MMio = diskmargin(P,C,sigma);

%% Results

fprintf('\n=== MULTILOOP DISK MARGINS ===\n');

fprintf('Input margin:        %.6f\n', MMi.DiskMargin);
fprintf('Output margin:       %.6f\n', MMo.DiskMargin);
fprintf('Concurrent I/O:      %.6f\n', MMio.DiskMargin);

fprintf('\n=== CONCURRENT I/O ===\n');
disp(MMio);

fprintf('\nNormalized comparison: 2*alpha_IO = %.6f\n', 2*MMio.DiskMargin);


% Gain-phase uncertainty region
figure('Name','LQI Gain-Phase Margin');
diskmarginplot(MMio.GainMargin);
grid on;

% Complex uncertainty disk
figure('Name','LQI Uncertainty Disk');
diskmarginplot(MMio.GainMargin,'disk');
grid on;

%% Disk Margins vs Frequency

% Input-only multiloop margins
Li = C*P;

% Output-only multiloop margins
Lo = P*C;

figure('Name','LQI Disk Margins vs Frequency');
diskmarginplot(Li, 'b', Lo, 'r');
legend('Input Margins', 'Output Margins');
grid on;

%% Pole Calculations

CL = Atilde - Btilde*KlqI;
poles = eig(CL);
[V,D] = eig(CL);

poles
%plot(real(poles), imag(poles),'rx')


tolerance = 1e-5; % Tune this


names = {'x','y','z','vx','vy','vz','P','Q','R', ...
         'p1','p2','p3','ex','ey','ez'};

for k = 1:15
    lambda = D(k,k);

    % Normalize eigenvector magnitudes for comparison
    participation = abs(V(:,k));
    participation = participation / max(participation);

    [sortedValues, idx] = sort(participation, 'descend');

    if real(lambda) < -1e-8
        stability = 'Stable';
    elseif real(lambda) > 1e-8
        stability = 'UNSTABLE';
    else
        stability = 'Marginal';
    end

    fprintf('\nMode %2d: lambda = %+.5f %+.5fi (%s)\n', ...
        k, real(lambda), imag(lambda), stability);

    if abs(imag(lambda)) > 1e-8
        fprintf('  Frequency: %.4f Hz, Period: %.4f s\n', ...
            abs(imag(lambda))/(2*pi), ...
            2*pi/abs(imag(lambda)));
    end

    fprintf('  Dominant states:\n');
    for j = 1:5
        fprintf('    %-3s  %.3f\n', ...
            names{idx(j)}, sortedValues(j));
    end
end


% Print eigenvalue and eigenvector phase information
for k = 1:15
    lambda = D(k,k);
    vnormalized = V(:,k);

    idx = find(abs(vnormalized) > tolerance);

    % Sort components by magnitude, largest first
    [~, order] = sort(abs(vnormalized(idx)), 'descend');
    idx = idx(order);

    phi = rad2deg(angle(vnormalized));
    parts = compose('%s(%+.0f°)', string(names(idx)).', phi(idx));

    fprintf('pole %8.3f%+8.3fi, phase=: %s\n', ...
        real(lambda), imag(lambda), strjoin(parts, ', '));
end

% Analyze and plot each eigenmode
tspan = linspace(0, 20, 1000);

for k = 1:15
    lambda = D(k,k);

    sigma = -real(lambda);
    wd = abs(imag(lambda));
    wn = abs(lambda);

    % Damping ratio
    if wn > 1e-12
        zeta = sigma / wn;
    else
        zeta = NaN;
    end

    % Oscillation period is undefined for a real eigenvalue
    if wd > 1e-12
        period = 2*pi/wd;
    else
        period = Inf;
    end

    % Settling and half-life times apply to stable modes
    if sigma > 1e-12
        ts = 4/sigma;
        thalf = log(2)/sigma;
    else
        ts = NaN;
        thalf = NaN;
    end

    % Exact real-valued response of this eigenmode
    % This also handles complex eigenvalues without ode45.
    X = zeros(numel(tspan), 15);
    for j = 1:numel(tspan)
        X(j,:) = real((V(:,k) * exp(lambda*tspan(j))).');
    end

    figure('Name', sprintf('Mode %d', k));

    subplot(4,1,1);
    plot(tspan, X(:,1:3));
    ylabel('pos [m]');
    legend('x','y','z');
    grid on;

    subplot(4,1,2);
    plot(tspan, X(:,4:6));
    ylabel('vel [m/s]');
    legend('vx','vy','vz');
    grid on;

    subplot(4,1,3);
    plot(tspan, X(:,7:9));
    ylabel('rates [rad/s]');
    legend('P','Q','R');
    grid on;

    subplot(4,1,4);
    plot(tspan, X(:,10:12));
    ylabel('attitude (MRP)');
    legend('p1','p2','p3');
    xlabel('t [s]');
    grid on;

    line1 = sprintf('Mode %d: \\lambda = %.3f %+.3fi', ...
        k, real(lambda), imag(lambda));

    line2 = sprintf( ...
        '\\zeta = %.3f | \\omega_n = %.4f rad/s | \\omega_d = %.4f rad/s | T = %.2f s', ...
        zeta, wn, wd, period);

    line3 = sprintf('t_s = %.2f s | t_{1/2} = %.2f s', ts, thalf);

    sgtitle({line1, line2, line3}, 'FontSize', 10);
end
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


%% Nyquist and Disk Margin Analysis of LQI Feedback

K = KlqI;

% Open-loop return ratio:
% L(s) = K * (sI - Atilde)^(-1) * Btilde
%
% Negative feedback:
% Acl = Atilde - Btilde*K

L = ss(Atilde, Btilde, K, ...
    zeros(size(K,1), size(Btilde,2)));

n = size(L,1);

assert(size(L,2) == n, 'Return ratio must be square.');

%% 1. Conventional MIMO Nyquist plot

figure('Name', 'Conventional MIMO Nyquist', ...
    'NumberTitle', 'off');

nyquist(L);
grid on;
sgtitle('LQI Return Ratio - Individual I/O Channels');


%% 2. Generalized MIMO Nyquist eigenloci
%% Generalized MIMO Nyquist visualization

% Use a more moderate initial frequency range
w = logspace(-3, 4, 5000);

Ljw = freqresp(L, w);
n = size(L,1);
N = numel(w);

eigL = zeros(n,N);
eigL(:,1) = eig(Ljw(:,:,1));

% Track eigenvalues by continuity
assignments = perms(1:n);

for k = 2:N
    currentEig = eig(Ljw(:,:,k));
    previousEig = eigL(:,k-1);

    bestCost = inf;

    for p = 1:size(assignments,1)
        candidate = currentEig(assignments(p,:));
        cost = sum(abs(candidate-previousEig).^2);

        if cost < bestCost
            bestCost = cost;
            bestEig = candidate;
        end
    end

    eigL(:,k) = bestEig;
end

% Diagnostics
fprintf('Maximum eigenvalue magnitude: %.3e\n', ...
    max(abs(eigL(:))));

fprintf('Minimum sampled distance from -1: %.3e\n', ...
    min(abs(eigL(:)+1)));

colors = lines(n);

figure('Name','Generalized MIMO Nyquist');

tiledlayout(1,2);

% Plot 1: Full frequency response
nexttile;
hold on; grid on;

for i = 1:n
    plot(real(eigL(i,:)), imag(eigL(i,:)), ...
        'Color',colors(i,:));
    plot(real(eigL(i,:)), -imag(eigL(i,:)), ...
        '--','Color',colors(i,:));
end

plot(-1,0,'rx','MarkerSize',12,'LineWidth',2);
title('Full Eigenvalue Response');
xlabel('Real');
ylabel('Imaginary');
axis equal;

% Plot 2: Critical region
nexttile;
hold on; grid on;

for i = 1:n
    plot(real(eigL(i,:)), imag(eigL(i,:)), ...
        'Color',colors(i,:), 'LineWidth',1.3);
    plot(real(eigL(i,:)), -imag(eigL(i,:)), ...
        '--','Color',colors(i,:));
end

plot(-1,0,'rx','MarkerSize',12,'LineWidth',2);

xlim([-5 3]);
ylim([-4 4]);
axis equal;

title('Zoom Near Critical Point');
xlabel('Real');
ylabel('Imaginary');



%% MIMO Disk Margin and Complex Gain Stability Map

K = KlqI;
A = Atilde;
B = Btilde;

% Open-loop return ratio
L = ss(A, B, K, zeros(size(K,1), size(B,2)));

%% 1. Calculate actual simultaneous MIMO disk margin

[DM, MM] = diskmargin(L);

alpha = MM.DiskMargin;
sigma = 0;  % Default balanced disk margin

disp('Simultaneous MIMO Disk Margin:');
disp(MM);

fprintf('Disk margin alpha: %.6f\n', alpha);
fprintf('Gain margin: [%.4f, %.4f]\n', MM.GainMargin);
fprintf('Phase margin: [%.2f, %.2f] deg\n', MM.PhaseMargin);

%% 2. Calculate the complex-gain stability map

realVals = linspace(-1.5, 3, 201);
imagVals = linspace(-2, 2, 201);

[Re, Im] = meshgrid(realVals, imagVals);

maxRealPole = zeros(size(Re));

BK = B*K;

for idx = 1:numel(Re)

    F = Re(idx) + 1i*Im(idx);

    % Apply common complex gain perturbation
    Acl = A - F*BK;

    polesCL = eig(Acl);

    maxRealPole(idx) = max(real(polesCL));
end

tol = 1e-9;

stable = maxRealPole < -tol;
unstable = maxRealPole > tol;

%% 3. Generate the actual MIMO uncertainty disk

theta = linspace(0, 2*pi, 1000);
delta = exp(1i*theta);

% Balanced disk from the actual multiloop margin
diskBoundary = ...
    (1 + alpha/2*delta) ./ ...
    (1 - alpha/2*delta);

%% 4. Plot the stability map and uncertainty disk

figure('Name','MIMO LQI Stability Region', ...
    'NumberTitle','off','Color','w');

hold on;
grid on;
box on;

% Stable points
scatter(Re(stable), Im(stable), 8, ...
    [0.1 0.65 0.25], 'filled', ...
    'DisplayName','Stable');

% Unstable points
scatter(Re(unstable), Im(unstable), 8, ...
    [0.85 0.1 0.15], 'filled', ...
    'DisplayName','Unstable');

% Marginal stability boundary
contour(Re, Im, maxRealPole, [0 0], ...
    'LineColor','m', ...
    'LineWidth',2, ...
    'HandleVisibility','off');

% Actual simultaneous MIMO disk
plot(real(diskBoundary), imag(diskBoundary), ...
    'k-', 'LineWidth',2.5, ...
    'DisplayName','MIMO Disk Margin');

% Nominal system
plot(1,0,'kx','MarkerSize',10, ...
    'LineWidth',2, ...
    'DisplayName','Nominal System');

xlabel('Real(F)');
ylabel('Imag(F)');

title(sprintf( ...
    'LQI Complex Gain Stability | MIMO Disk Margin = %.4f', ...
    alpha));

axis equal;
xlim([min(realVals),max(realVals)]);
ylim([min(imagVals),max(imagVals)]);

legend('Location','best');
hold off;

%% 5. MATLAB native disk margin visualizations

% Gain-phase tradeoff
figure('Name','MIMO Gain-Phase Margin');
diskmarginplot(MM.GainMargin);

% Actual uncertainty disk
figure('Name','MIMO Uncertainty Disk');
diskmarginplot(MM.GainMargin,'disk');

%% Figure 3: Margins vs Frequency
figure('Name','Disk Margins vs Frequency');

diskmarginplot(L);
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
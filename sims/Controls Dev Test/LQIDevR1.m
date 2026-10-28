% clear all;



% states = [x  3]

%xmaxs = [1 1 1 1 1 1 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4];
%umaxs = [20 pi/16 pi/16 5];

%Qdiags = 1./xmaxs.^2;
%Rdiags = 1./umaxs.^2;


Qdiags = [0 0 0 10 10 10 1e2 1e2 1e2 1e6 1e6 1e6];
QIdiags =[10 10 500];
QI = diag(horzcat(Qdiags,QIdiags));
Q = diag(Qdiags);

% T delta rcs
Rdiags =[100,1e8,1e8,1e4];

R = diag(Rdiags);


params = struct();
params.g = 9.81;

params.d  = .2;
params.r_m = 0.125;
params.n_rcs = 1;
params.I = [1 0 0; 0 1 0; 0 0 1];

m = 1;

q = [sqrt(2)/2;0;sqrt(2)/2;0];
%q=[0 0 0 0]';
p = qToMRP(q)

x = [0;0;0;0;0;0;0;0;0;p];
u  = [10;0;0;0];

[A, B] = HopperLinearization_lqi(x, u, params, m);



%%

[Atilde,Btilde] = extendedDynamics(A,B);
CrossTerm = zeros(15,4);

%sysI = sys(Atilde,Btilde,CI,DI);

KlqI = lqr(Atilde,Btilde,QI,R,CrossTerm)
% KlqI(3,:) = KlqI(2,:);
% KlqI(3,2)= KlqI(2,1);
% KlqI(3,5)= KlqI(2,4);
% KlqI(3,9) = KlqI(2,8);
% KlqI(3,10) = KlqI(2,11);
% KlqI(3,14) = KlqI(2,13);
%KlqI(3,12) = -KlqI(2,11);%Klqi = lqi(sys,Q,R,CrossTerm);
% %KlqI(:,10:12)

%KlqI(:,10:12)


%% Nyquist stability analysis of the LQI feedback loop

% Open-loop return ratio:
% L(s) = K * (sI - Atilde)^(-1) * Btilde
%
% Negative feedback gives:
% Acl = Atilde - Btilde*K

K = KlqI;

% Construct the 4-input, 4-output return-ratio system.
% Its output is K*x_aug, and its input is the control input.
L = ss(Atilde, Btilde, K, zeros(size(K,1), size(Btilde,2)));

%% 1. Conventional Nyquist plot
figure('Name', 'LQI Nyquist Plot');
nyquist(L);
grid on;
title('Nyquist Plot of LQI Return Ratio');



%% Generalized MIMO Nyquist eigenloci

w = logspace(-4, 2, 3000);
Ljw = freqresp(L, w);

eigL = zeros(size(L,1), length(w));

for k = 1:length(w)
    eigL(:,k) = eig(Ljw(:,:,k));
end

% Explicitly create a new figure and axes
figNyq = figure( ...
    'Name', 'Generalized MIMO Nyquist Plot', ...
    'NumberTitle', 'off');

axNyq = axes('Parent', figNyq);
hold(axNyq, 'on');
grid(axNyq, 'on');

% Plot each eigenvalue locus
for k = 1:size(eigL,1)
    plot(axNyq, real(eigL(k,:)), imag(eigL(k,:)), ...
        'LineWidth', 1.2);
    plot(axNyq, real(eigL(k,:)), -imag(eigL(k,:)), ...
        '--', 'HandleVisibility', 'off');
end

% Critical point
plot(axNyq, -1, 0, 'rx', ...
    'MarkerSize', 12, 'LineWidth', 2);

xlabel(axNyq, 'Real');
ylabel(axNyq, 'Imaginary');
title(axNyq, 'Generalized MIMO Nyquist Eigenloci');
legend(axNyq, {'Eigenlocus 1', 'Eigenlocus 2', ...
    'Eigenlocus 3', 'Eigenlocus 4', 'Critical point -1'}, ...
    'Location', 'best');

hold(axNyq, 'off');
drawnow;



%% 3. MIMO disk margins (if supported by your MATLAB version)
try
    DM = diskmargin(L);
    disp('MIMO disk-margin analysis:');
    disp(DM);
    figure('Name', 'Disk Margin');
    diskmargin(L);
catch ME
    warning('Disk-margin analysis unavailable: %s', ME.message);
end




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

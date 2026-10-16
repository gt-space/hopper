 clear all; close all; clc;


% states = [x y z vx vy vz P Q R p1 p2 p3]

%xmaxs = [1 1 1 1 1 1 1e-3 1e-3 1e-3 1e-4 1e-4 1e-4];
%umaxs = [20 pi/16 pi/16 5];

%Qdiags = 1./xmaxs.^2;
%Rdiags = 1./umaxs.^2;


Qdiags = [0 0 0 30 30 30 1 1 1 10 10 10];
QIdiags =[10 10 500];
QI = diag(horzcat(Qdiags,QIdiags));
Q = diag(Qdiags);

% T delta rcs
Rdiags =[100,1e9,1e9,1];

R = diag(Rdiags);


params = struct();
params.g = 9.81;

params.d  = .2;
params.r_m = 0.125;
params.n_rcs = 1;
params.I = [1 0 0; 0 4 0; 0 0 4];

m = 1;

q = [sqrt(2)/2;0;sqrt(2)/2;0];
%q=[1 0 0 0]';
p = qToMRP(q)

x = [0;0;0;0;0;0;0;0;0;p];
u  = [10;1e-4;1e-4;0];

[A, B] = HopperLinearization_lqi(x, u, params, m);



%%

[Atilde,Btilde] = extendedDynamics(A,B);
CrossTerm = zeros(15,4);

%sysI = sys(Atilde,Btilde,CI,DI);

KlqI = lqr(Atilde,Btilde,QI,R,CrossTerm)
%Klqi = lqi(sys,Q,R,CrossTerm);
% KlqI(:,10:12)



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
  vnormalized = V(:,k);
  idx    = find(vnormalized > tolerance);
  [vals, order] = sort(vnormalized(idx), 'descend');
  idx    = idx(order);
  %parts = arrayfun(@(j) sprintf('%s = %.3f', names{j}, v(j)), idx, 'UniformOutput', false);
  phi = rad2deg(angle(V(:,k)));
  parts = compose('%s(%+.0f°)', string(names(idx)).', phi(idx));
  fprintf('pole %8.3f%+8.3fi, phase= : %s\n', real(D(k,k)), imag(D(k,k)), strjoin(parts, ', '))
end
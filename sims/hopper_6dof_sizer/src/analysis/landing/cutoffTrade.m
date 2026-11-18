function X = cutoffTrade(mcFile, opts)
%CUTOFFTRADE  Engine cutoff height x crush-pad plateau: P(pads don't crush).
%   X = cutoffTrade(hopperRoot('results','mc_touchdown_200.mat'))
%   X = cutoffTrade(X)          re-evaluate a saved trade (no new sims)
%
%   1. From each flight Monte Carlo run (runMonteCarlo), the state when the
%      engine cut (descent rate v_c, thrust) and at touchdown (drift, tilt,
%      rates, mass; touchdownCases).
%   2. For each run, pad plateau and surface, the touchdown sink rate at
%      which the pads start to crush (bisection on landingSim, engine off),
%      and whether they crush when the engine is still on at contact.
%   3. For each commanded cutoff height h and altitude-knowledge error
%      sigma: the actual cut height is h + e, e ~ N(0, sigma). Above the
%      ground the vehicle free-falls (the flight model cuts thrust
%      instantly): sink = sqrt(v_c^2 + 2 g (h + e)). At or below it, the
%      feet touch first with the engine on at the descent rate v_c, and
%      touchdown detection shuts it down: thrust decays from contact with
%      time constant tau (ShutdownTau, the trade's third axis).
%   P(no crush) is over runs and e samples. Drift and tilt are taken from
%   the touchdown state at the flight model's own 5 cm cutoff.
%
%   X.P(h, sigma, tau, plateau, surface)  P(no crush)
%   X.Pon(tau, plateau, surface)          P(no crush), engine on at contact
%   X.Preq(h, sigma)                 P(sink <= the 1.5 m/s requirement)
%   X.vstar(run, plateau, surface)   crush-onset sink rate (m/s)
%   Saves results/cutoff_trade.mat and outputs/landing/cutoff_trade.png.
%   Each crush-onset search is appended to LogFile as it finishes (with a
%   progress line); rerunning skips searches already in the log.
%
%   With ShutdownProfile, thrust after the cutoff command follows that curve
%   (fraction of the thrust at the command) instead of dropping to zero:
%   the fall from the cut height is integrated with it, the landing starts
%   with the thrust left at contact, and the search is on the cut height
%   (X.vstar then holds the highest cut height that does not crush).
%   Engine-on landings cut Latency after contact; the tau axis becomes the
%   latency axis.

arguments
    mcFile
    opts.Plateau double = [4 5 6 8] * 1e3
    opts.Surfaces string = ["concrete", "soil_firm"]
    opts.Heights double = -0.06:0.005:0.10            % m, commanded cutoff height (< 0: only touchdown detection cuts, unless the height error is larger)
    opts.Sigma double = [0 0.01 0.02 0.03]            % m, cutoff-height error (1 sigma)
    opts.ShutdownTau double = [0.05 0.2]              % s, thrust decay after contact (engine on at contact)
    opts.PlotTau (1,1) double = 0.05                  % s, tau shown in the figure
    opts.Bisect (1,1) double = 4                      % bisection steps on 0.2-1.6 m/s (0.09 m/s bracket)
    opts.SimTime (1,1) double = 0.4                   % s simulated per landing
    opts.LogFile (1,:) char = fullfile(hopperRoot(), 'results', 'cutoff_trade_log.csv')
    opts.FromLog (1,1) logical = false                % only read the log (interim result: runs with every search done)
    opts.ShutdownProfile double = []                  % [t; fraction] thrust after the cutoff command (e.g. from the
                                                      % fluid network); replaces the instant cut and ShutdownTau
    opts.Latency double = [0 0.05]                    % s, touchdown detection to cutoff command (with a profile)
    opts.CrushTol (1,1) double = 0.5e-3               % m of pad stroke that counts as crushed
    opts.Samples (1,1) double = 400                   % e samples per run
    opts.Parallel (1,1) logical = true
    opts.OutFile (1,:) char = fullfile(hopperRoot(), 'results', 'cutoff_trade.mat')
    opts.FigDir (1,:) char = fullfile(hopperRoot(), 'outputs', 'landing')
end

if isstruct(mcFile)                % saved trade: re-evaluate only
    X = mcFile;
else
    X = thresholds(mcFile, opts);
end
X = evaluate(X, opts);
save(opts.OutFile, 'X');
plotTrade(X, opts.FigDir);
end

% ---------------------------------------------------------------------
function X = thresholds(mcFile, opts)
L = load(mcFile);
VC = vehicleConfig('hopper');
[cases, info] = touchdownCases(L.results, Vehicle=VC);
ok = find(info.used);
res = L.results(ok);
nR = numel(cases);
vc = nan(1, nR);  Tc = nan(1, nR);  hc = nan(1, nR);
for k = 1:nR
    ts = res(k).timeseries;
    i = find(ts.thrust > 0.5 * max(ts.thrust(ts.t > ts.t(end) - 1)), 1, 'last');
    vc(k) = ts.vz_raw(i);  Tc(k) = ts.thrust(i);  hc(k) = ts.altitude(i);
end
fprintf('cutoffTrade: %d flight runs; descent at cutoff %.3f +- %.3f m/s, cut at %.3f m\n', ...
    nR, mean(vc), std(vc), median(hc));

nF = numel(opts.Plateau);  nS = numel(opts.Surfaces);  nT = numel(opts.ShutdownTau);
[ff, ss, kk] = ndgrid(1:nF, 1:nS, 1:nR);      % every configuration for run 1, then run 2, ...
nJob = numel(kk);
plat = opts.Plateau;  surfs = opts.Surfaces;  tau = opts.ShutdownTau;  tol = opts.CrushTol;
nb = opts.Bisect;  tEnd = opts.SimTime;
prof = opts.ShutdownProfile;
if ~isempty(prof)
    tau = opts.Latency;  nT = numel(tau);       % engine-on axis is detection latency
end

% Each finished search is appended to the log at once, so a run can be
% watched and resumed; jobs already in the log are skipped.
done = readLog(opts.LogFile, nT);
todo = setdiff(1:nJob, done.job);
if opts.FromLog, todo = []; end
if isempty(done.job)
    fid = fopen(opts.LogFile, 'w');
    fprintf(fid, 'job,run,plateau_kN,surface,vstar%s\n', sprintf(',crush_on_tau%g', tau));
    fclose(fid);
end
fprintf('cutoffTrade: %d crush-onset searches, %d already logged; log: %s\n', nJob, nJob - numel(todo), opts.LogFile);
t0 = tic;
q = parallel.pool.DataQueue;
afterEach(q, @(r) logRow(r, opts.LogFile, numel(todo), t0, plat, surfs));
logRow([], '', 0, 0, [], []);                 % reset the counter
if isempty(todo)
    % nothing to run (FromLog, or all logged)
elseif opts.Parallel
    parfor jj = 1:numel(todo)
        j = todo(jj);
        [v, on] = onset(cases(kk(j)), vc(kk(j)), Tc(kk(j)), plat(ff(j)), surfs(ss(j)), tau, tol, VC, nb, tEnd, prof);
        send(q, [j, kk(j), ff(j), ss(j), v, on(:).']);
    end
else
    for jj = 1:numel(todo)
        j = todo(jj);
        [v, on] = onset(cases(kk(j)), vc(kk(j)), Tc(kk(j)), plat(ff(j)), surfs(ss(j)), tau, tol, VC, nb, tEnd, prof);
        logRow([j, kk(j), ff(j), ss(j), v, on(:).'], opts.LogFile, numel(todo), t0, plat, surfs);
    end
end
fprintf('cutoffTrade: searches done in %.0f s\n', toc(t0));

D = readLog(opts.LogFile, nT);
fi = arrayfun(@(F) find(plat == F * 1e3, 1), D.plateau);
si = arrayfun(@(s) find(surfs == s, 1), D.surface);
X.vstar = nan(nR, nF, nS);
X.complete = false(1, nR);
X.crush_engine_on = false(nR, nT, nF, nS);
for r = 1:numel(D.job)
    X.vstar(D.run(r), fi(r), si(r)) = D.vstar(r);
    X.crush_engine_on(D.run(r), :, fi(r), si(r)) = D.on(r, :);
end
full = accumarray(D.run(:), 1, [nR 1]).' == nF * nS;
X.complete = full;
X.vstar = X.vstar(full, :, :);
X.crush_engine_on = X.crush_engine_on(full, :, :, :);
vc = vc(full);  Tc = Tc(full);  hc = hc(full);
fprintf('cutoffTrade: %d of %d runs have every search done\n', nnz(full), nR);
X.tau = tau;
X.profile = prof;
X.vc = vc;  X.Tc = Tc;  X.hc = hc;
X.info = structfun(@(a) a(ok(full)), rmfield(info, 'used'), 'uni', 0);
X.plateau = plat;  X.surfaces = surfs;
X.mcFile = mcFile;
end

function logRow(r, file, nTodo, t0, plat, surfs)
persistent n
if isempty(r), n = 0; return; end
n = n + 1;
fid = fopen(file, 'a');
fprintf(fid, '%d,%d,%g,%s,%.4f%s\n', r(1), r(2), plat(r(3)) / 1e3, surfs(r(4)), r(5), sprintf(',%d', r(6:end)));
fclose(fid);
el = toc(t0);
fprintf('[%4d/%d] run %3d  %4.1f kN  %-9s  crush onset %5.2f m/s  engine-on crush %s   %5.1f min, ~%4.1f min left\n', ...
    n, nTodo, r(2), plat(r(3)) / 1e3, surfs(r(4)), r(5), mat2str(r(6:end)), el / 60, el / 60 / n * (nTodo - n));
end

function D = readLog(file, nT)
D = struct('job', [], 'run', [], 'plateau', [], 'surface', strings(0, 1), 'vstar', [], 'on', zeros(0, nT));
if ~isfile(file), return; end
T = readtable(file, 'TextType', 'string', 'VariableNamingRule', 'modify');
if isempty(T), return; end
T = T(all(~ismissing(T(:, [1:3 5:5 + nT])), 2) & ~ismissing(T.surface), :);   % drop a row cut short by an interrupted write
D.job = T.job;  D.run = T.run;  D.plateau = T.plateau_kN;  D.surface = T.surface;  D.vstar = T.vstar;
D.on = logical(T{:, 6:5 + nT});
end

function [vs, onCrush] = onset(c, vc, Tc, F, surface, tau, tol, VC, nb, tEnd, prof)
if ~isempty(prof)
    [vs, onCrush] = onsetProfile(c, vc, Tc, F, surface, tau, tol, VC, nb, tEnd, prof);
    return
end
% Crush-onset sink rate with the engine off (bisection; the largest rate
% seen not to crush, so conservative by up to one bracket), and whether
% the pads crush with the engine on at contact for each shutdown tau.
G = landingGear(Vehicle=VC, Core=crushCore('pad', Plateau=F));
S = landingSurface(char(surface));
crushed = @(v) run(G, S, setfield(setfield(c, 'sink', v), 'thrust', 0), tol, tEnd); %#ok<SFLD>
lo = 0.2;  hi = 1.6;                 % 1.6 m/s is above any sink in the trade (10 cm cut: 1.44)
if crushed(lo), vs = 0;
elseif ~crushed(hi), vs = Inf;
else
    for it = 1:nb
        mid = (lo + hi) / 2;
        if crushed(mid), hi = mid; else, lo = mid; end
    end
    vs = lo;
end
onCrush = false(numel(tau), 1);
for t = 1:numel(tau)
    ce = c;  ce.sink = vc;  ce.thrust = Tc;  ce.thrust_tau = tau(t);
    onCrush(t) = run(G, S, ce, tol, tEnd);
end
end

function [hs, onCrush] = onsetProfile(c, vc, Tc, F, surface, lat, tol, VC, nb, tEnd, prof)
% Highest cut height (m) whose landing does not crush (bisection on
% 0-12 cm; the largest height seen not to crush), and engine-on crush for
% each detection latency.
G = landingGear(Vehicle=VC, Core=crushCore('pad', Plateau=F));
S = landingSurface(char(surface));
c.thrust = Tc;  c.thrust_profile = prof;
crushedAt = @(h) cutLanding(G, S, c, vc, h, tol, tEnd, prof);
lo = 0;  hi = 0.12;
if crushedAt(lo), hs = -Inf;              % crushes even when cut at contact
elseif ~crushedAt(hi), hs = Inf;
else
    for it = 1:nb
        mid = (lo + hi) / 2;
        if crushedAt(mid), hi = mid; else, lo = mid; end
    end
    hs = lo;
end
onCrush = false(numel(lat), 1);
for t = 1:numel(lat)
    ce = c;  ce.sink = vc;  ce.thrust_t0 = -lat(t);
    onCrush(t) = run(G, S, ce, tol, tEnd);
end
end

function tf = cutLanding(G, S, c, vc, h, tol, tEnd, prof)
[tFall, v] = fall(h, vc, c.thrust, c.mass, prof);
c.sink = v;  c.thrust_t0 = tFall;
tf = run(G, S, c, tol, tEnd);
end

function [t, v] = fall(h, vc, T0, m, prof)
% Fall from height h at descent rate vc with thrust T0 * prof(t) (up).
g = 9.80665;  dt = 1e-4;  t = 0;  z = h;  v = vc;          % v positive down
while z > 0 && t < 5
    f = interp1(prof(1, :), prof(2, :), min(max(t, prof(1, 1)), prof(1, end)));
    a = g - T0 * f / m;
    v = v + a * dt;  z = z - v * dt;  t = t + dt;
end
v = max(v, 0);
end

function tf = run(G, S, c, tol, tEnd)
R = landingSim(G, S, c, Log=false, tEnd=tEnd);
tf = R.peak.stroke_pad > tol;
end

% ---------------------------------------------------------------------
function X = evaluate(X, opts)
g = 9.80665;
rng(7);
e = randn(1, opts.Samples);
nH = numel(opts.Heights);  nG = numel(opts.Sigma);
if ~isfield(X, 'profile'), X.profile = []; end
nR = numel(X.vc);  nF = numel(X.plateau);  nS = numel(X.surfaces);  nT = numel(X.tau);
P = zeros(nH, nG, nT, nF, nS);  Preq = zeros(nH, nG);  v95 = zeros(nH, nG);
for a = 1:nH
    for b = 1:nG
        h = opts.Heights(a) + opts.Sigma(b) * e;              % 1 x M
        air = h > 0;
        v = sqrt(X.vc(:).^2 + 2 * g * max(h, 0));             % nR x M
        Preq(a, b) = mean(v(:) <= 1.5);
        v95(a, b) = prctile(v(:), 95);
        for t = 1:nT
            for f = 1:nF
                for s = 1:nS
                    if isempty(X.profile)
                        okAir = v < X.vstar(:, f, s);
                    else
                        okAir = h < X.vstar(:, f, s);            % vstar holds the crush-onset cut height
                    end
                    ok = okAir .* air + (~X.crush_engine_on(:, t, f, s)) .* ~air;
                    P(a, b, t, f, s) = mean(ok(:));
                end
            end
        end
    end
end
X.heights = opts.Heights;  X.sigma = opts.Sigma;
X.P = P;  X.Preq = Preq;  X.v95 = v95;
X.Pon = squeeze(mean(~X.crush_engine_on, 1));            % tau x plateau x surface
if isempty(X.profile)
    X.Pair = squeeze(mean(X.vc(:) < X.vstar, 1));        % P(no crush) if the cut came at the instant of contact, engine off
end
X.plot_tau = opts.PlotTau;
end

% ---------------------------------------------------------------------
function plotTrade(X, figDir)
[~, it] = min(abs(X.tau - X.plot_tau));
if ~exist(figDir, 'dir'), mkdir(figDir); end
f = figure('Visible', 'off', 'Position', [0 0 1500 850], 'Color', 'w');
try, theme(f, 'light'); catch, end
nG = numel(X.sigma);
tl = tiledlayout(f, 1, nG, 'TileSpacing', 'compact');
if ~isempty(X.profile)
    title(tl, sprintf(['Runs where the pads do not crush vs commanded cutoff height, fluid-network shutdown ' ...
        '(solid: concrete, dashed: firm soil; engine-on landings cut %g s after contact)'], X.tau(it)));
else
title(tl, sprintf(['Runs where the pads do not crush vs commanded engine cutoff height ' ...
    '(solid: concrete, dashed: firm soil; engine-on landings shut down with \\tau = %g s)'], X.tau(it)));
end
co = lines(numel(X.plateau));
for b = 1:nG
    ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
    for fi = 1:numel(X.plateau)
        for s = 1:numel(X.surfaces)
            ls = '-';  if s > 1, ls = '--'; end
            plot(ax, X.heights * 100, 100 * X.P(:, b, it, fi, s), ls, 'Color', co(fi, :), 'LineWidth', 1.4, ...
                'DisplayName', sprintf('%g kN %s', X.plateau(fi) / 1e3, strrep(X.surfaces(s), '_', ' ')), ...
                'HandleVisibility', ternary(s == 1, 'on', 'off'));
        end
    end
    yline(ax, 95, 'k:', '95 %', 'HandleVisibility', 'off');
    xlabel(ax, 'commanded cutoff height (cm)');
    if b == 1, ylabel(ax, 'runs with no crush (%)'); end
    title(ax, sprintf('cutoff-height error \\sigma = %g cm', 100 * X.sigma(b)));
    ylim(ax, [0 100]);
    if b == nG, legend(ax, 'Location', 'southwest'); end
end
exportgraphics(f, fullfile(figDir, 'cutoff_trade.png'), 'Resolution', 150);
close(f);
end

function o = ternary(c, a, b)
if c, o = a; else, o = b; end
end

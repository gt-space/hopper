function C = tvcLoopCampaign(n, opts)
%TVCLOOPCAMPAIGN  Model-in-the-loop TVC trade: Monte Carlo per actuator config.
%   C = tvcLoopCampaign(50) flies the nominal hopper plus n seeded Monte
%   Carlo scenarios (mc_params.json, the same scenarios for every config)
%   in models/hopper_6dof.slx for each TVC configuration below, with the
%   TVC_MODEL variant, in parallel (parsim). Saves results/tvc_loop.mat.
%
%   Per run: landing and requirement results (mc_main + checkConstraints,
%   as runMonteCarlo scores them) and the gimbal behaviour: commanded
%   (u_cmd) and achieved (act_bus) angle, rate, acceleration, tracking
%   error, and the share of flight time spent at the rate limit.
%
%   Ordinary dispersions barely load the gimbal, so with Stress=true
%   (default) every config also flies three deterministic cases built on
%   the nominal scenario: a 12 m/s crosswind, a 2.5 cm thrust-line offset
%   on both axes (the vehicle needs gimbal trim and the actuators carry the
%   offset moment), and both together. n = 0 flies the nominal (and the
%   stress cases) only.
%
%   A config may set the drive's clamps, as fields of opts.Configs:
%     accel_dps2  setpoint acceleration ramp (deg/s^2), [] or Inf = off
%     rate_dps    setpoint speed clamp (deg/s),         [] or Inf = off
%   (tvcWorkspace AccelLimit / RateLimit). Per run the row then also gives
%   peak achieved rate/accel over the clamp, the share of flight time at
%   each clamp, attitude error (tilt from vertical) and a sustained
%   gimbal-oscillation check (see oscillation()).

arguments
    n (1,1) double = 50
    opts.Configs struct = defaultConfigs()
    opts.OutputFile (1,:) char = hopperRoot('results', 'tvc_loop.mat')
    opts.Workers (1,1) double = NaN
    opts.Stress (1,1) logical = true
end
nom = nominalScenario();
scen = {nom};  names = "nominal";
if n > 0                                      % n = 0: nominal (+ stress) only
    scen = [scen, num2cell(generateScenarios('mc_params.json', n))];
    names = [names, "mc_" + string(1:n)];
end
if opts.Stress
    w = 12 / sqrt(2);                         % 12 m/s resultant, diagonal
    s1 = nom;  s1.uwind = w;  s1.vwind = w;
    s2 = nom;  s2.engine_off_axis_y = 0.025;  s2.engine_off_axis_z = 0.025;
    s3 = s2;   s3.uwind = w;  s3.vwind = w;
    scen = [scen, {s1, s2, s3}];
    names = [names, "stress_wind12", "stress_offset25", "stress_wind12_offset25"];
end
cfg = opts.Configs;
runs = {};
for c = 1:numel(cfg)
    for k = 1:numel(scen)
        runs{end+1} = struct('cfg', cfg(c), 'c', c, 'k', k, 'case', names(k), 'scenario', scen{k}); %#ok<AGROW>
    end
end
fprintf('TVC loop campaign: %d configs x %d scenarios = %d runs\n', numel(cfg), numel(scen), numel(runs));

if isempty(gcp('nocreate'))
    nw = opts.Workers;  if isnan(nw), nw = mcWorkerCount(numel(runs)); end
    parpool('Processes', nw);
end
ins = repmat(Simulink.SimulationInput('hopper_6dof'), 1, numel(runs));
for j = 1:numel(runs)
    r = runs{j};
    ins(j) = ins(j).setPreSimFcn(@(in) prepare(in, r)) ...
                   .setPostSimFcn(@(out) summarize(out, r));
end
t0 = tic;
so = parsim(ins, 'ShowProgress', 'on', 'TransferBaseWorkspaceVariables', 'off', ...
            'UseFastRestart', 'off', 'StopOnError', 'off');
C.elapsed = toc(t0);
C.configs = cfg;
C.n = n;
R = repmat(emptyRow(), numel(runs), 1);
for j = 1:numel(runs)
    o = so(j);
    row = [];
    try, row = o.get('row'); catch, end
    if isempty(row)
        row = emptyRow();
        row.config = string(cfg(runs{j}.c).name);  row.scenario = runs{j}.k;  row.case = runs{j}.case;
        row.error = string(o.ErrorMessage);
    end
    R(j) = row;
end
% Runs whose model failed to compile (parallel workers occasionally collide
% building the same MATLAB Function blocks) are retried once, serially.
retry = find(arrayfun(@(r) ~r.success && contains(r.error, ["compile", "propagat", "determine sizes"]), R));
for j = retry'
    fprintf('retrying run %d (%s, scenario %d) serially\n', j, R(j).config, R(j).scenario);
    try
        o = sim(prepare(Simulink.SimulationInput('hopper_6dof'), runs{j}));
        R(j) = getfield(summarize(o, runs{j}), 'row'); %#ok<GFLD>
    catch e
        R(j).error = string(e.message);
    end
end
C.retried = numel(retry);
C.runs = struct2table(R);
C.summary = groupsummary(C.runs(startsWith(C.runs.case, ["nominal", "mc_"]), :), 'config', {'mean', 'max'}, ...
    {'landing_offset', 'max_tilt', 'max_cmd_deg', 'max_rate_dps', 'max_accel_dps2', 'rms_err_deg', 'max_err_deg', 'rate_sat_pct'});
save(opts.OutputFile, '-struct', 'C');
fprintf('saved %s (%.0f s)\n', opts.OutputFile, C.elapsed);
end

%% --------------------------------------------------------------------
function cfg = defaultConfigs()
cfg = struct('name', {}, 'model', {}, 'layout', {}, 'offset', {});
cfg(end+1) = struct('name', "PDR 6DOF (2nd order, 15 deg/s)", 'model', "second_order", 'layout', "", 'offset', []);
cfg(end+1) = struct('name', "AG 9113 recommended", 'model', "ag9113", 'layout', "recommended", 'offset', []);
cfg(end+1) = struct('name', "AG 9313 recommended", 'model', "ag9313", 'layout', "recommended", 'offset', []);
cfg(end+1) = struct('name', "AG 9313 sizer layout", 'model', "ag9313", 'layout', "sizer", 'offset', []);
cfg(end+1) = struct('name', "AG 9113 recommended, 2.5 cm load", 'model', "ag9113", 'layout', "recommended", 'offset', 0.025);
cfg(end+1) = struct('name', "AG 9313 sizer layout, 2.5 cm load", 'model', "ag9313", 'layout', "sizer", 'offset', 0.025);
end

function in = prepare(in, r)
V = variantConfig('TVC_MODEL', r.cfg.model, 'LOG', fullfile(tempdir, sprintf('tvcloop_%d_%d.mat', r.c, r.k)));
if isempty(r.scenario)
    W = hopperWorkspace([], V);
else
    W = hopperWorkspace(r.scenario, V);
end
if r.cfg.model ~= "second_order"
    [alim, rlim] = clamps(r.cfg);
    W = tvcWorkspace(W, V, Layout=char(r.cfg.layout), Offset=r.cfg.offset, ...
                     AccelLimit=deg2rad(alim), RateLimit=deg2rad(rlim));
end
in = buildSimInput(in, [], W, V);
end

function S = summarize(out, r)
S = struct();
row = emptyRow();
row.config = string(r.cfg.name);  row.scenario = r.k;  row.case = r.case;
logFile = fullfile(tempdir, sprintf('tvcloop_%d_%d.mat', r.c, r.k));
cleanup = onCleanup(@() deleteIf(logFile));
try
    if ~isempty(out.ErrorMessage), error('run:failed', '%s', out.ErrorMessage); end
    % ---- landing and requirements, as runMonteCarlo scores them ----
    res = mc_main(r.scenario, out);
    row.success = res.status.success;
    if res.status.success
        req = checkConstraints(res);
        row.meets_requirements = req.status.pass;
        row.missed = strjoin(string(req.status.failed), ", ");
        row.landing_offset = req.flight.landing_offset;
        row.landing_vel = req.flight.landing_vel;
        row.max_tilt = req.flight.max_tilt_deg;
        row.touchdown_tilt = req.flight.touchdown_tilt_deg;
    end
    % ---- gimbal ----
    L = load(logFile, 'logsout');
    act = L.logsout.get('act_bus').Values;
    t = act.tvc_pitch.Time;
    ach = [squeeze(act.tvc_pitch.Data), squeeze(act.tvc_yaw.Data)];
    u = L.logsout.get('gnc_bus').Values.u_cmd;     % [thrust, pitch, yaw, rcs]
    D = squeeze(u.Data);
    if size(D, 1) ~= numel(u.Time), D = D.'; end   % time down the rows
    cmd = interp1(u.Time, D(:, 2:3), t, 'previous', 'extrap');
    [t, iu] = unique(t);  ach = ach(iu, :);  cmd = cmd(iu, :);
    tt = (t(1):1e-3:t(end))';
    ach = interp1(t, ach, tt);  cmd = interp1(t, cmd, tt, 'previous');
    [b, a] = butter(2, 25 / 500);                 % 25 Hz: keeps the actuator's own dynamics
    af = filtfilt(b, a, ach);
    rate = gradient(af', 1e-3)';  acc = gradient(rate', 1e-3)';
    err = cmd - ach;
    row.max_cmd_deg = rad2deg(max(abs(cmd(:))));
    row.max_ach_deg = rad2deg(max(abs(ach(:))));
    row.max_rate_dps = rad2deg(max(abs(rate(:))));
    row.p99_rate_dps = rad2deg(prctile(abs(rate(:)), 99));
    row.max_accel_dps2 = rad2deg(max(abs(acc(:))));
    row.rms_err_deg = rad2deg(sqrt(mean(err(:).^2)));
    row.max_err_deg = rad2deg(max(abs(err(:))));
    % time at the rate limit: |rate| within 90 % of the largest rate the
    % configuration reached anywhere in the run, while error > 0.2 deg
    lim = 0.9 * max(abs(rate(:)));
    row.rate_sat_pct = 100 * mean(any(abs(rate) >= lim & abs(err) > deg2rad(0.2), 2));
    % ---- against the drive clamps (achieved gimbal, 25 Hz filtered) ----
    % at clamp: within 5 % of it on either axis
    [alim, rlim] = clamps(r.cfg);
    row.accel_lim_dps2 = alim;  row.rate_lim_dps = rlim;
    if isfinite(rlim)
        row.rate_over_lim = row.max_rate_dps / rlim;
        row.at_rate_lim_pct = 100 * mean(any(abs(rad2deg(rate)) >= 0.95 * rlim, 2));
    end
    if isfinite(alim)
        row.accel_over_lim = row.max_accel_dps2 / alim;
        row.at_accel_lim_pct = 100 * mean(any(abs(rad2deg(acc)) >= 0.95 * alim, 2));
    end
    % ---- attitude error: tilt from vertical (the hop's reference) ----
    if res.status.success
        tilt = abs(90 - res.timeseries.pitch);
        row.rms_tilt_deg = sqrt(mean(tilt.^2, 'omitnan'));
    end
    % ---- sustained gimbal oscillation ----
    O = oscillation(tt, ach);
    row.osc_amp_deg = O.amp_deg;  row.osc_freq_hz = O.freq_hz;
    row.osc_dur_s = O.dur_s;      row.osc_sustained = O.sustained;
catch e
    row.success = false;
    row.error = string(e.message);
end
S.row = row;
end

function row = emptyRow()
row = struct('config', "", 'scenario', NaN, 'case', "", 'success', false, 'meets_requirements', false, ...
    'missed', "", 'landing_offset', NaN, 'landing_vel', NaN, 'max_tilt', NaN, 'touchdown_tilt', NaN, ...
    'max_cmd_deg', NaN, 'max_ach_deg', NaN, 'max_rate_dps', NaN, 'p99_rate_dps', NaN, ...
    'max_accel_dps2', NaN, 'rms_err_deg', NaN, 'max_err_deg', NaN, 'rate_sat_pct', NaN, ...
    'rate_lim_dps', Inf, 'accel_lim_dps2', Inf, 'rate_over_lim', NaN, 'accel_over_lim', NaN, ...
    'at_rate_lim_pct', NaN, 'at_accel_lim_pct', NaN, 'rms_tilt_deg', NaN, ...
    'osc_amp_deg', NaN, 'osc_freq_hz', NaN, 'osc_dur_s', NaN, 'osc_sustained', false, 'error', "");
end

function [alim, rlim] = clamps(cfg)
% drive clamps of a config, deg/s^2 and deg/s (Inf = off)
alim = Inf;  rlim = Inf;
if isfield(cfg, 'accel_dps2') && ~isempty(cfg.accel_dps2), alim = cfg.accel_dps2; end
if isfield(cfg, 'rate_dps') && ~isempty(cfg.rate_dps), rlim = cfg.rate_dps; end
end

function O = oscillation(t, ach)
% Sustained gimbal oscillation: achieved angle band-passed 0.5-10 Hz, its
% amplitude (sqrt(2) x RMS) over 1 s windows stepped 0.25 s. Sustained =
% the amplitude stays above 0.1 deg on one axis for 3 s or more in a row.
% The first and last 1 s (lift-off, touchdown) are left out. Reports the
% longest such stretch (dur_s), the largest window amplitude, and the
% dominant frequency of that stretch (or of the whole flight if none).
thr = 0.1;  minDur = 3;  win = 1;  hop = 0.25;  fs = 1 / (t(2) - t(1));
[b, a] = butter(2, [0.5 10] / (fs / 2), 'bandpass');
x = rad2deg(filtfilt(b, a, ach));
k = t >= t(1) + 1 & t <= t(end) - 1;  t = t(k);  x = x(k, :);
starts = t(1):hop:(t(end) - win);
amp = zeros(numel(starts), size(x, 2));
for j = 1:numel(starts)
    m = t >= starts(j) & t < starts(j) + win;
    amp(j, :) = sqrt(2) * sqrt(mean(x(m, :).^2, 1));
end
O.amp_deg = max(amp(:));
best = 0;  ax = 1;  i0 = 1;
for c = 1:size(amp, 2)
    on = [false; amp(:, c) > thr; false];
    d = diff(on);  s = find(d == 1);  e = find(d == -1) - 1;
    [L, i] = max(e - s + 1);
    if ~isempty(L) && L > best, best = L;  ax = c;  i0 = s(i); end
end
O.dur_s = 0;  m = true(size(t));
if best > 0
    O.dur_s = (best - 1) * hop + win;
    m = t >= starts(i0) & t < starts(i0) + O.dur_s;
end
O.sustained = O.dur_s >= minDur;
X = abs(fft(x(m, ax) - mean(x(m, ax))));
f = (0:numel(X) - 1)' * fs / numel(X);
h = f > 0.5 & f < 10;
[~, i] = max(X(h));  fh = f(h);
O.freq_hz = NaN;  if ~isempty(i), O.freq_hz = fh(i); end
end

function s = nominalScenario()
p = readstruct('mc_params.json');
s = struct();
for f = fieldnames(p)'
    s.(f{1}) = p.(f{1}).nominal;
end
end

function deleteIf(f)
if isfile(f), delete(f); end
end

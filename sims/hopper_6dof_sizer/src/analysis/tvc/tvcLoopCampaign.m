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
%   offset moment), and both together.

arguments
    n (1,1) double = 50
    opts.Configs struct = defaultConfigs()
    opts.OutputFile (1,:) char = hopperRoot('results', 'tvc_loop.mat')
    opts.Workers (1,1) double = NaN
    opts.Stress (1,1) logical = true
end
nom = nominalScenario();
scen = [{nom}, num2cell(generateScenarios('mc_params.json', n))];
names = ["nominal", "mc_" + string(1:n)];
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
    alim = Inf;
    if isfield(r.cfg, 'accel_dps2') && ~isempty(r.cfg.accel_dps2), alim = deg2rad(r.cfg.accel_dps2); end
    W = tvcWorkspace(W, V, Layout=char(r.cfg.layout), Offset=r.cfg.offset, AccelLimit=alim);
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
    'max_accel_dps2', NaN, 'rms_err_deg', NaN, 'max_err_deg', NaN, 'rate_sat_pct', NaN, 'error', "");
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

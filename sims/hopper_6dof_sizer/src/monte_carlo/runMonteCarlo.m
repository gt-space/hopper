function R = runMonteCarlo(n, opts)
%RUNMONTECARLO  Seeded Monte Carlo of a hopper model, serial or with parsim.
%   R = runMonteCarlo(10)
%   R = runMonteCarlo(100, UseParallel=true, OutputFile=hopperRoot('results','mc_100.mat'))
%
%   Options
%     Model        model name (default hopper_6dof_NED_v2)
%     ParamsFile   dispersion bounds (default mc_params.json)
%     Variables    struct of extra model variables, e.g. struct('enable_slosh', 0)
%     UseParallel  run with parsim. Uses the open pool if there is one;
%                  otherwise starts one sized to fit in free RAM (each
%                  worker commits ~3.3 GB; see mcWorkerCount)
%     Workers      pool size to start instead of the automatic choice
%     Nominal      also run the nominal scenario (default true)
%     Limits       overrides for checkConstraints
%     OutputFile   save R here ("" = don't save)
%
%   R.results   1xn mc_main results with status.pass / status.failed (against
%               Limits) and status.meets_requirements / missed_requirements
%               (against the mission requirements in inputs/mission_inputs.m)
%   R.nominal   nominal result (or [])
%   R.summary   table, one row per run
%   R.scenarios the generated scenarios (generateScenarios, rng seed 42)

arguments
    n (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.Model (1,:) char = 'hopper_6dof_NED_v2'
    opts.ParamsFile (1,:) char = 'mc_params.json'
    opts.Variables (1,1) struct = struct()
    opts.UseParallel (1,1) logical = false
    opts.Workers (1,1) double = NaN
    opts.Nominal (1,1) logical = true
    opts.Limits (1,1) struct = struct()
    opts.OutputFile (1,1) string = ""
end

scenarios = generateScenarios(opts.ParamsFile, n);
runs = num2cell(scenarios);
if opts.Nominal
    runs = [{nominalScenario(opts.ParamsFile)}, runs];
end

t0 = tic;
fprintf('Monte Carlo: %d scenarios of %s (%s)\n', n, opts.Model, ...
    ternary(opts.UseParallel, 'parsim', 'serial'));
if opts.UseParallel
    outs = runParallel(runs, opts);
else
    outs = cell(size(runs));
    for k = 1:numel(runs)
        fprintf('  run %d / %d\n', k, numel(runs));
        try
            outs{k} = sim(prepareInput(Simulink.SimulationInput(opts.Model), runs{k}, opts.Variables));
        catch err
            outs{k} = err;
        end
    end
end

results = cell(size(runs));
for k = 1:numel(runs)
    results{k} = postProcess(runs{k}, outs{k}, opts.Limits);
end
if opts.Nominal
    R.nominal = results{1};
    results = results(2:end);
else
    R.nominal = [];
end
R.results = toStructArray(results);
R.scenarios = scenarios;
R.summary = summaryTable(R.results);
R.elapsed = toc(t0);
R.model = opts.Model;

fprintf('\n%s\n', repmat('=', 1, 60));
disp(R.summary);
fprintf('%d / %d passed, %d errored, %d met mission requirements, %.0f s\n', nnz(R.summary.pass), n, ...
    nnz(~R.summary.success), nnz(R.summary.meets_requirements), R.elapsed);
if strlength(opts.OutputFile) > 0
    save(opts.OutputFile, '-struct', 'R');
    fprintf('saved %s\n', opts.OutputFile);
end
end

function in = prepareInput(in, scenario, extra)
W = hopperWorkspace(scenario);
for f = fieldnames(extra)'
    W.(f{1}) = extra.(f{1});
end
in = buildSimInput(in, [], W);
end

function outs = runParallel(runs, opts)
% parsim opens the current project on each worker, so the path is there.
if isempty(gcp('nocreate'))
    nw = opts.Workers;
    if isnan(nw), nw = mcWorkerCount(numel(runs)); end
    fprintf('Starting a parallel pool with %d workers\n', nw);
    parpool('Processes', nw);
end
ins = repmat(Simulink.SimulationInput(opts.Model), size(runs));
for k = 1:numel(runs)
    s = runs{k};
    ins(k) = ins(k).setPreSimFcn(@(in) prepareInput(in, s, opts.Variables));
end
so = parsim(ins, 'ShowProgress', 'on', 'TransferBaseWorkspaceVariables', 'off');
outs = num2cell(so);
end

function r = postProcess(scenario, out, limits)
if isa(out, 'MException')
    r = emptyResult(scenario, out.message);
    return
end
if ~isempty(out.ErrorMessage)
    r = emptyResult(scenario, out.ErrorMessage);
    return
end
r = mc_main(scenario, out);
if r.status.success
    req = checkConstraints(r);              % mission requirements
    r = checkConstraints(r, limits);        % limits this run is judged on
    r.status.meets_requirements = req.status.pass;
    r.status.missed_requirements = req.status.failed;
else
    r.status.pass = false;
    r.status.failed = {'diverged'};
    r.status.meets_requirements = false;
    r.status.missed_requirements = {'diverged'};
end
end

function r = emptyResult(scenario, msg)
r.scenario = scenario;
r.status.success = false;
r.status.error = msg;
r.status.pass = false;
r.status.failed = {'error'};
r.status.meets_requirements = false;
r.status.missed_requirements = {'error'};
end

function S = toStructArray(c)
% Errored runs carry fewer fields; pad them so the runs concatenate.
f = cellfun(@fieldnames, c, 'uni', 0);
names = unique(vertcat(f{:}), 'stable');
for k = 1:numel(c)
    for f = setdiff(names, fieldnames(c{k}))'
        c{k}.(f{1}) = [];
    end
    c{k} = orderfields(c{k}, names);
end
S = [c{:}];
end

function T = summaryTable(results)
n = numel(results);
run = (1:n)';
[success, pass, meets_requirements] = deal(false(n, 1));
[max_alt, landing_vel, landing_offset, max_tilt, flight_time] = deal(nan(n, 1));
[failed, missed_requirements] = deal(strings(n, 1));
for k = 1:n
    r = results(k);
    success(k) = r.status.success;
    pass(k) = r.status.pass;
    meets_requirements(k) = r.status.meets_requirements;
    failed(k) = strjoin(string(r.status.failed), ", ");
    missed_requirements(k) = strjoin(string(r.status.missed_requirements), ", ");
    if success(k)
        max_alt(k) = r.flight.max_altitude;
        landing_vel(k) = r.flight.landing_vel;
        landing_offset(k) = r.flight.landing_offset;
        max_tilt(k) = r.flight.max_tilt_deg;
        flight_time(k) = r.flight.flight_time;
    else
        failed(k) = string(r.status.error);
    end
end
T = table(run, success, pass, max_alt, landing_vel, landing_offset, max_tilt, flight_time, failed, ...
    meets_requirements, missed_requirements);
end

function s = nominalScenario(paramsFile)
p = readstruct(paramsFile);
s = struct();
for f = fieldnames(p)'
    s.(f{1}) = p.(f{1}).nominal;
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

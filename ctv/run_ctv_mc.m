% RUN_CTV_MC  Run the CTV Monte Carlo end to end. Open this file and press Run.
%
% Works from a fresh MATLAB in any folder: it moves to ctv/ itself, prints the
% parameters it is about to use, runs the batch, plots, and summarises.
%
% ---------------------------------------------------------------------------
% Editing the parameters: mc_params_ctv.json is the only place to set nominal
% values, bounds and distributions. One entry per parameter, for example:
%
%   "ox_mass":      {"lower": 15, "upper": 17, "nominal": 16}
%   "mass_factor":  {"type": "normal", "mean": 1.0, "std": 0.01, "nominal": 1.0,
%                    "lower": 0.95, "upper": 1.05}
%   "cstar_eff":    {"type": "lognormal", "mean": 0, "std": 0.05, "nominal": 1}
%   "gust":         {"type": "weibull", "scale": 5, "shape": 2, "nominal": 4}
%   "thrust_bias":  {"type": "beta", "alpha": 2, "beta": 5,
%                    "lower": -0.02, "upper": 0.02, "nominal": 0}
%   "vent_setting": {"type": "discrete", "values": [1, 2, 3], "nominal": 2}
%
% "type" defaults to uniform when absent. "nominal" is required on every entry
% and is what the nominal run uses. For every type except beta, "lower" and
% "upper" are optional and clip the samples.
%
% To PIN a parameter, set lower = upper (the nominal run still uses "nominal").
%
% To ADD a parameter, add an entry here AND a line in the Monte Carlo block at
% the bottom of inputs/mission_inputs.m that copies it onto the matching IN.*
% field, otherwise it is sampled but never reaches the model.
%
% To REMOVE a parameter, delete its entry; the model then uses whatever
% mission_inputs.m sets for it by default.
%
% Changing n or the parameter file invalidates the checkpoints, which are
% cleared automatically.
% ---------------------------------------------------------------------------

% ======================== Settings ========================
n          = 10;                   % number of Monte Carlo runs
runMode    = 'serial';             % 'serial' | 'serial_isolated' | 'nominal'
makePlots  = true;
paramsFile = 'mc_params_ctv.json';
% ==========================================================

cd(fileparts(mfilename('fullpath')));
addpath(pwd);

show_params(paramsFile);

mc_parallelization_ctv   % reads n, runMode, makePlots, paramsFile

report_run();


function show_params(paramsFile)
%SHOW_PARAMS Print the parameter set the run is about to use.
if isempty(fileparts(paramsFile))
    paramsFile = fullfile(pwd, paramsFile);
end
bounds = readstruct(paramsFile);
names = fieldnames(bounds);

fprintf('\nParameters from %s\n', paramsFile);
fprintf('%-26s %-10s %10s  %-30s %s\n', 'name', 'dist', 'nominal', 'spread', 'pinned');
fprintf('%s\n', repmat('-', 1, 92));
for k = 1:numel(names)
    p = bounds.(names{k});
    if isfield(p, 'type'); ty = char(p.type); else; ty = 'uniform'; end

    if isfield(p, 'nominal'); nomStr = sprintf('%10g', p.nominal);
    else;                     nomStr = sprintf('%10s', '(missing)');
    end

    switch ty
        case 'uniform'
            spread = sprintf('[%g, %g]', p.lower, p.upper);
        case {'normal', 'lognormal'}
            spread = sprintf('mean %g, std %g', p.mean, p.std);
        case 'weibull'
            spread = sprintf('scale %g, shape %g', p.scale, p.shape);
        case 'beta'
            spread = sprintf('a %g, b %g on [%g, %g]', p.alpha, p.beta, p.lower, p.upper);
        case 'discrete'
            spread = ['{' strjoin(arrayfun(@(v) sprintf('%g', v), p.values(:).', ...
                'UniformOutput', false), ', ') '}'];
        otherwise
            spread = '(unrecognised type)';
    end

    pinned = 'no';
    if all(isfield(p, {'lower', 'upper'})) && isequal(p.lower, p.upper)
        pinned = 'yes';
    end

    fprintf('%-26s %-10s %s  %-30s %s\n', names{k}, ty, nomStr, spread, pinned);
end
fprintf('\n');
end


function report_run()
%REPORT_RUN Summarise the newest results file in this folder.
d = dir(fullfile(pwd, 'mc_results_ctv_*.mat'));
if isempty(d)
    fprintf('No results file was produced.\n');
    return
end
[~, newest] = max([d.datenum]);
resultsFile = fullfile(d(newest).folder, d(newest).name);

S = load(resultsFile);
if ~isfield(S, 'results')
    fprintf('Nominal-only run.\nResults file: %s\n', resultsFile);
    return
end

R = S.results;
nr = numel(R);
ok    = arrayfun(@(r) r.status.success, R);
pass  = arrayfun(@(r) r.status.pass, R);
passF = arrayfun(@(r) r.status.pass_flight, R);

fprintf('%s\n', repmat('=', 1, 62));
fprintf('Runs: %d | errors: %d | pass: %d | pass_flight: %d\n', ...
    nr, sum(~ok), sum(pass), sum(passF));

if any(ok)
    good  = R(ok);
    drift = arrayfun(@(r) max(hypot(r.timeseries.x, r.timeseries.y)), good);
    rows = { ...
        'landing velocity (m/s)', arrayfun(@(r) r.flight.landing_vel,  good); ...
        'max altitude (m)',       arrayfun(@(r) r.flight.max_altitude, good); ...
        'drift (m)',              drift};
    fprintf('%s\n', repmat('-', 1, 62));
    fprintf('%-24s %10s %10s %10s\n', 'metric', 'mean', 'p5', 'p95');
    for k = 1:size(rows, 1)
        v = rows{k, 2};
        fprintf('%-24s %10.3f %10.3f %10.3f\n', ...
            rows{k, 1}, mean(v), prctile(v, 5), prctile(v, 95));
    end
end
fprintf('%s\n', repmat('=', 1, 62));
fprintf('Results file:   %s\n', resultsFile);
fprintf('Scenario table: %s\n', fullfile(pwd, 'mc_scenarios.csv'));
end

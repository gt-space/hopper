function results = mc_runner(n_scenarios, params_file, output_file, model, flags)
% mc_runner(n_scenarios, params_file, output_file, model, flags)
%   Serial Monte Carlo, kept for its old call signature; runMonteCarlo does
%   the work (and takes more options, e.g. UseParallel).
%   model: the old speed copies, now harness configurations of hopper_6dof:
%     'original'   (default) every axis at its default
%     'toggleable' slosh and wind switched by flags
%     'minimal'    SLOSH=off, MASS_PROPS=constant, no wind
%   flags: [enable_slosh enable_wind] for 'toggleable' (1 = on), default [1 1]
%   Pass [] to keep any default, e.g.  mc_runner(20, [], [], 'minimal')

if nargin < 1 || isempty(n_scenarios); n_scenarios = 10;              end
if nargin < 2 || isempty(params_file); params_file = 'mc_params.json'; end
if nargin < 4 || isempty(model);       model       = 'original';       end
if nargin < 5 || isempty(flags);       flags       = [1 1];            end
if nargin < 3 || isempty(output_file)
    output_file = hopperRoot('results', 'mc_results.mat');
    if ~strcmpi(model, 'original'); output_file = hopperRoot('results', sprintf('mc_results_%s.mat', lower(model))); end
end

switch lower(model)
    case 'original',   V = variantConfig(); wind = true;
    case 'toggleable'
        V = variantConfig();
        if ~flags(1), V = variantConfig('SLOSH', 'off'); end
        wind = logical(flags(2));
    case 'minimal',    V = variantConfig('SLOSH', 'off', 'MASS_PROPS', 'constant'); wind = false;
    otherwise
        error('mc_runner:unknownModel', ...
            'Unknown model "%s". Use ''original'', ''toggleable'' or ''minimal''.', model);
end
if ~strcmpi(model, 'toggleable') && ~isequal(flags(:)', [1 1])
    warning('mc_runner:flagsIgnored', 'Flags only apply to ''toggleable''; ignoring them for "%s".', model);
end

R = runMonteCarlo(n_scenarios, Variants=V, Wind=wind, ParamsFile=params_file, OutputFile=output_file);
results = R.results;
end

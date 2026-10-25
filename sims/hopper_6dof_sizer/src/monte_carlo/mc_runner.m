function results = mc_runner(n_scenarios, params_file, output_file, model, flags)
% mc_runner(n_scenarios, params_file, output_file, model, flags)
%   Serial Monte Carlo, kept for its old call signature; runMonteCarlo does
%   the work (and takes more options, e.g. UseParallel).
%   model: 'original' (default) | 'toggleable' | 'minimal'   (see mc_sim_input.m)
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

% mc_sim_input picks the model and sets its toggle variables
in = mc_sim_input(model, flags);
extra = struct();
for v = in.Variables
    extra.(v.Name) = v.Value;
end

R = runMonteCarlo(n_scenarios, Model=in.ModelName, ParamsFile=params_file, ...
    Variables=extra, OutputFile=output_file);
results = R.results;
end

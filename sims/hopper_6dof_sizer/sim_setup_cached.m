% sim_setup_cached  Fast workspace setup for the 6DOF models.
%
% First run (or after sim_setup.m / inputs change): runs the full sim_setup
% and saves the resulting workspace to sim_workspace_cache.mat.
% Later runs: loads that cache in seconds instead of re-running everything.
%
% Force a rebuild with:   rebuild_cache = true; sim_setup_cached
%
% Also defines the enable_slosh / enable_wind toggles used by
% hopper_6dof_NED_v2_simplified (1 = feature on, 0 = feature off).
% They default to 1 (full physics). Turn a feature off before simulating, e.g.
%   enable_slosh.Value = 0;
% or, for a single run without touching the workspace,
%   in = Simulink.SimulationInput('hopper_6dof_NED_v2_simplified');
%   out = sim(in.setVariable('enable_slosh', 0));
%
% hopper_6dof_NED_v2_minimal computes its own frozen mass / inertia / CG values
% in its InitFcn, so it needs nothing extra from this script.

here       = fileparts(mfilename('fullpath'));
cache_file = fullfile(here, 'sim_workspace_cache.mat');
cd(here);

if exist('rebuild_cache', 'var') && rebuild_cache
    if isfile(cache_file), delete(cache_file); end
end
clear rebuild_cache

if isfile(cache_file)
    fprintf('Loading cached workspace: %s\n', cache_file);
    load(cache_file);
else
    fprintf('No cache found, running full sim_setup (slow, one time)...\n');
    tic;
    sim_setup;                       % note: sim_setup begins with clearvars
    here       = fileparts(which('sim_setup'));
    cache_file = fullfile(here, 'sim_workspace_cache.mat');
    fprintf('sim_setup finished in %.1f s, saving cache...\n', toc);
    save(cache_file, '-v7.3');
end

% Feature toggles (1 = feature on, 0 = off)
enable_slosh  = Simulink.Parameter(1);
enable_wind   = Simulink.Parameter(1);
for p = {enable_slosh, enable_wind}
    p{1}.DataType = 'double';
    p{1}.StorageClass = 'Auto';
end
clear p

fprintf('Workspace ready.\n');

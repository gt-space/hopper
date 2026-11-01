% sim_setup_cached  Fast base-workspace setup (used by fsw_bridge).
%
% First run (or after sim_setup.m / inputs change): runs the full sim_setup
% and saves the resulting workspace to sim_workspace_cache.mat.
% Later runs: loads that cache in seconds instead of re-running everything.
%
% Force a rebuild with:   rebuild_cache = true; sim_setup_cached

cache_file = hopperRoot('work', 'sim_workspace_cache.mat');
cd(hopperRoot);

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
    cache_file = hopperRoot('work', 'sim_workspace_cache.mat');
    fprintf('sim_setup finished in %.1f s, saving cache...\n', toc);
    save(cache_file, '-v7.3');
end

fprintf('Workspace ready.\n');

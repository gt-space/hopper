function fixBernoulliSeeds()
%FIXBERNOULLISEEDS  Give every wind-gust Bernoulli generator its own seed.
%   The gust on/off Bernoulli Binary Generator copied into the 6DOF models
%   had SeedSource = 'Auto', which ignores its Seed and draws from MATLAB's
%   global random stream. Results then depended on whatever used the RNG
%   earlier in the session (a previous sim, generateScenarios' rng(42), test
%   order), so two identical runs in one session differed. This sets
%   SeedSource = 'Parameter', keeping the Seed already on the block, in every
%   model under models/, and saves the ones it changed.

f = dir(hopperRoot('models', '**', '*.slx'));
for k = 1:numel(f)
    [~, m] = fileparts(f(k).name);
    wasLoaded = bdIsLoaded(m);
    load_system(fullfile(f(k).folder, f(k).name));
    b = find_system(m, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'IncludeCommented', 'on', ...
        'MatchFilter', @Simulink.match.allVariants, 'MaskType', 'Bernoulli Binary Generator');
    b = b(strcmp(get_param(b, 'SeedSource'), 'Auto'));
    for j = 1:numel(b)
        set_param(b{j}, 'SeedSource', 'Parameter');
        fprintf('%s: %s now uses Seed = %s\n', m, erase(b{j}, [m '/']), get_param(b{j}, 'Seed'));
    end
    if ~isempty(b), save_system(m); end
    if ~wasLoaded, close_system(m, 0); end
end
end

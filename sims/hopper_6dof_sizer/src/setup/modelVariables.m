function names = modelVariables(model, W)
%MODELVARIABLES  Names of the base-workspace variables a model reads.
%   names = modelVariables(model, W) compiles the model with the workspace
%   struct W (from hopperWorkspace) and asks Simulink.findVars which
%   variables it uses. The answer is cached per model until the .slx changes.
%
%   For a model with variants (models/hopper_6dof.slx) only the active
%   choices are compiled, so the answer depends on the VARIANT_* values in
%   W and is cached per combination. The VARIANT_* control variables
%   themselves are always included.

persistent cache
if isempty(cache), cache = containers.Map(); end

file = which(model);
if isempty(file), error('modelVariables:noModel', 'Model %s is not on the path.', model); end
info = dir(file);
vf = fieldnames(W);
vf = sort(vf(startsWith(vf, 'VARIANT_')));
key = sprintf('%s|%s%s', model, info.date, strjoin(cellfun(@(f) sprintf('|%s=%g', f, W.(f)), vf, 'uni', 0), ''));
if isKey(cache, key)
    names = cache(key);
    return
end

% findVars compiles the model, so put W in the base workspace for the
% duration (the caller's base variables are restored afterwards).
saved = evalin('base', 'who');
savedVals = cellfun(@(n) evalin('base', n), saved, 'uni', 0);
restore = onCleanup(@() restoreBase(saved, savedVals));
evalin('base', 'clearvars');
for f = fieldnames(W)'
    assignin('base', f{1}, W.(f{1}));
end

wasLoaded = bdIsLoaded(model);
load_system(model);
if ~wasLoaded, closeModel = onCleanup(@() close_system(model, 0)); end
vars = Simulink.findVars(model, 'SourceType', 'base workspace', 'SearchMethod', 'compiled');
names = unique([{vars.Name}, variantControls(model)]);
cache(key) = names;
end

function names = variantControls(model)
% VARIANT_* variables named in any variant choice's control expression
blks = find_system(model, 'MatchFilter', @Simulink.match.allVariants, 'LookUnderMasks', 'all', ...
    'FollowLinks', 'on', 'RegExp', 'on', 'VariantControl', 'VARIANT_');
exprs = get_param(blks, 'VariantControl');
names = unique([regexp(strjoin(cellstr(exprs), ' '), 'VARIANT_\w+', 'match'), {}]);
end

function restoreBase(names, vals)
evalin('base', 'clearvars');
for k = 1:numel(names)
    assignin('base', names{k}, vals{k});
end
end

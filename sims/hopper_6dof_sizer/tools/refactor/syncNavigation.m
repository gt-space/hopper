function syncNavigation(branch)
%SYNCNAVIGATION  Merge the nav team's latest filter work into this branch.
%   syncNavigation()            merge origin/GMEKF
%   syncNavigation('origin/x')  merge another branch of theirs
%
%   The nav team works under sims/hopper_6dof_sizer/NavigationSims/ (the
%   layout before the MATLAB project); here their files live in
%   subsystems/navigation/. This fetches, merges with --no-commit, moves
%   anything that landed under NavigationSims/ into subsystems/navigation/
%   and drops files they deleted.
%
%   Some of their files are edited on this branch too (the requests in
%   docs/VARIANTS.md, Navigation: diagnostics out of the correct step, the
%   torque process noise, the fixed-step CTV filter, nav_sensor_params...).
%   A file changed on this branch since the last merge of theirs is never
%   resolved automatically: if both sides changed it, it is left
%   conflicted for a person. Every other conflict takes their side, and
%   those files are checked to match their branch exactly.
%
%   It stops, merge left open, if the merge touches anything outside the
%   navigation folders or leaves conflicts to resolve. Commit the merge
%   yourself after looking at it, then run syncProjectFiles.

if nargin < 1, branch = 'origin/GMEKF'; end
root = fullfile(hopperRoot(), '..', '..');   % repository root: git paths below are relative to it
git = @(cmd) runGit(root, cmd);
if ~isempty(strtrim(git('status --porcelain --untracked-files=no')))
    error('syncNavigation:dirty', 'Commit or set aside your changes first.');
end
git('fetch --quiet origin');
old = 'sims/hopper_6dof_sizer/NavigationSims/';
new = 'sims/hopper_6dof_sizer/subsystems/navigation/';

% files this branch changed since their last merged commit (ours to keep)
% (their commit has the files under the old folder: compare file by file)
base = strtrim(git(sprintf('merge-base HEAD %s', branch)));
here = splitlines(strtrim(git(sprintf('ls-tree -r --name-only HEAD "%s"', new))));
ours = {};
for k = 1:numel(here)
    if isempty(here{k}), continue; end
    a = strtrim(git(sprintf('rev-parse "HEAD:%s"', here{k})));
    b = strtrim(runGit(root, sprintf('rev-parse "%s:%s"', base, strrep(here{k}, new, old)), true));
    if ~strcmp(a, b), ours{end+1} = here{k}; end %#ok<AGROW>
end

[~] = runGit(root, sprintf('merge --no-ff --no-commit %s', branch), true);
st = splitlines(strtrim(git('status --porcelain')));
left = {};
for k = 1:numel(st)
    line = st{k};
    if strlength(line) < 4, continue; end
    code = line(1:2); path = strtrim(erase(line(4:end), '"'));
    if ~startsWith(path, {old, new, 'sims/hopper_6dof_sizer/resources/project/'}) && ~strcmp(path, '.gitignore')
        error('syncNavigation:outside', ['The merge changes %s, outside the navigation folders. ' ...
            'Resolve it by hand (the merge is still open).'], path);
    end
    if any(strcmp(code, {'UA', 'AU', 'AA', 'UU', 'UD', 'DU'})) && any(strcmp(path, ours))
        left{end+1} = path; %#ok<AGROW>
        continue
    end
    switch code
        case {'UA', 'AU', 'AA', 'UU'}, git(sprintf('checkout --theirs -- "%s"', path)); git(sprintf('add -- "%s"', path));
        case {'UD', 'DU'},             git(sprintf('rm -q -- "%s"', path));
    end
end
% anything still under the old folder moves to the new one
moved = splitlines(strtrim(git(sprintf('ls-files "%s"', old))));
for k = 1:numel(moved)
    if isempty(moved{k}), continue; end
    dst = strrep(moved{k}, old, new);
    d = fileparts(fullfile(root, dst));
    if ~isfolder(d), mkdir(d); end
    git(sprintf('mv -f "%s" "%s"', moved{k}, dst));
end
if ~isempty(left)
    error('syncNavigation:conflicts', ['Both sides changed these; resolve them by hand (the merge is ' ...
        'still open), keeping this branch''s changes where they still apply:\n  %s'], strjoin(left, '\n  '));
end

theirs = splitlines(strtrim(git(sprintf('ls-tree -r --name-only %s %s', branch, old))));
bad = {};
for k = 1:numel(theirs)
    f = theirs{k};
    here = strrep(f, old, new);
    if any(strcmp(here, ours)), continue; end   % edited here: merged, not copied
    a = strtrim(git(sprintf('rev-parse "%s:%s"', branch, f)));
    b = strtrim(runGit(root, sprintf('rev-parse ":%s"', here), true));
    if ~strcmp(a, b), bad{end+1} = here; end %#ok<AGROW>
end
if ~isempty(bad)
    error('syncNavigation:mismatch', 'These do not match %s: %s', branch, strjoin(bad, ', '));
end
fprintf('subsystems/navigation matches %s apart from this branch''s %d edited file(s). Review, commit the merge, then run syncProjectFiles.\n', ...
    branch, numel(ours));
end

function out = runGit(root, cmd, allowFail)
if nargin < 3, allowFail = false; end
[s, out] = system(sprintf('git -C "%s" %s', root, cmd));
if s ~= 0 && ~allowFail
    error('syncNavigation:git', 'git %s failed:\n%s', cmd, out);
end
end

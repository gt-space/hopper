function syncNavigation(branch)
%SYNCNAVIGATION  Merge the nav team's latest filter work into this branch.
%   syncNavigation()            merge origin/GMEKF
%   syncNavigation('origin/x')  merge another branch of theirs
%
%   The nav team works under sims/hopper_6dof_sizer/NavigationSims/ (the
%   layout before the MATLAB project); here their files live in
%   subsystems/navigation/. This fetches, merges with --no-commit, moves
%   anything that landed under NavigationSims/ into subsystems/navigation/,
%   takes their side for their files, drops files they deleted, and checks
%   that subsystems/navigation matches their branch file for file. It
%   stops, merge left open, if the merge touches anything outside the
%   navigation folders: that needs a person. Commit the merge yourself
%   after looking at it, then run syncProjectFiles.

if nargin < 1, branch = 'origin/GMEKF'; end
root = fullfile(hopperRoot(), '..', '..');   % repository root: git paths below are relative to it
git = @(cmd) runGit(root, cmd);
if ~isempty(strtrim(git('status --porcelain --untracked-files=no')))
    error('syncNavigation:dirty', 'Commit or set aside your changes first.');
end
git('fetch --quiet origin');
[~] = runGit(root, sprintf('merge --no-ff --no-commit %s', branch), true);

old = 'sims/hopper_6dof_sizer/NavigationSims/';
new = 'sims/hopper_6dof_sizer/subsystems/navigation/';
st = splitlines(strtrim(git('status --porcelain')));
for k = 1:numel(st)
    line = st{k};
    if strlength(line) < 4, continue; end
    code = line(1:2); path = strtrim(erase(line(4:end), '"'));
    if ~startsWith(path, {old, new, 'sims/hopper_6dof_sizer/resources/project/'}) && ~strcmp(path, '.gitignore')
        error('syncNavigation:outside', ['The merge changes %s, outside the navigation folders. ' ...
            'Resolve it by hand (the merge is still open).'], path);
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

theirs = splitlines(strtrim(git(sprintf('ls-tree -r --name-only %s %s', branch, old))));
bad = {};
for k = 1:numel(theirs)
    f = theirs{k};
    a = strtrim(git(sprintf('rev-parse "%s:%s"', branch, f)));
    b = strtrim(runGit(root, sprintf('rev-parse ":%s"', strrep(f, old, new)), true));
    if ~strcmp(a, b), bad{end+1} = strrep(f, old, new); end %#ok<AGROW>
end
if ~isempty(bad)
    error('syncNavigation:mismatch', 'These do not match %s: %s', branch, strjoin(bad, ', '));
end
fprintf('subsystems/navigation matches %s. Review, commit the merge, then run syncProjectFiles.\n', branch);
end

function out = runGit(root, cmd, allowFail)
if nargin < 3, allowFail = false; end
[s, out] = system(sprintf('git -C "%s" %s', root, cmd));
if s ~= 0 && ~allowFail
    error('syncNavigation:git', 'git %s failed:\n%s', cmd, out);
end
end

function plan = buildfile
%BUILDFILE  HopperSim build tasks (run from the project root with buildtool).
%   buildtool check      repo guards: no build caches / autosaves committed,
%                        no model saved newer than the team release
%                        (TEAM_RELEASE file)
%   buildtool test       unit, model-file and FDIR tests (fast)
%   buildtool sim        nominal flight and slosh open-loop tests
%   buildtool mc         seeded 10-run Monte Carlo against tests/montecarlo/mcLimits
%   buildtool            check + test
%   CI runs: buildtool check test sim mc
%
%   JUnit results go to outputs/test-results/ for CI.

import matlab.buildtool.tasks.TestTask

plan = buildplan(localfunctions);
results = fullfile('outputs', 'test-results');

plan("test") = TestTask(["tests/unit", "tests/models", "tests/fdir"], ...
    TestResults = fullfile(results, "test.xml"), ...
    Description = "Unit, model-file and FDIR tests");
plan("sim") = TestTask("tests/sim", ...
    TestResults = fullfile(results, "sim.xml"), ...
    Description = "Nominal flight and slosh open-loop tests");
plan("mc") = TestTask("tests/montecarlo", ...
    TestResults = fullfile(results, "mc.xml"), ...
    Description = "Seeded 10-run Monte Carlo");

plan("test").Dependencies = "check";
plan.DefaultTasks = ["check", "test"];
end

function checkTask(~)
% Repo guards that need no simulation.
root = hopperRoot();
[status, tracked] = system(sprintf('git -C "%s" ls-files', root));
assert(status == 0, 'git ls-files failed: %s', tracked);
tracked = splitlines(strtrim(tracked));

junk = tracked(~cellfun(@isempty, regexp(tracked, ...
    '(^|/)slprj/|\.slxc$|\.asv$|\.autosave$|\.slx\.r20\d\d[ab]$', 'once')));
assert(isempty(junk), 'Build caches / autosaves are committed:\n  %s', strjoin(junk, '\n  '));

team = teamRelease();
running = ['R' version('-release')];
if isNewerRelease(team, running)
    warning('buildfile:olderMATLAB', ...
        'This is MATLAB %s but the team release is %s; install %s.', running, team, team);
end
live = tracked(endsWith(tracked, '.slx') & ~startsWith(tracked, 'archive/'));
newer = {};
for k = 1:numel(live)
    rel = Simulink.MDLInfo(fullfile(root, live{k})).ReleaseName;
    if isNewerRelease(rel, team)
        newer{end+1} = sprintf('%s (%s)', live{k}, rel); %#ok<AGROW>
    end
end
assert(isempty(newer), 'Models saved newer than the team release %s:\n  %s', team, strjoin(newer, '\n  '));
fprintf('check: %d tracked files, %d live models, none newer than %s\n', numel(tracked), numel(live), team);
end

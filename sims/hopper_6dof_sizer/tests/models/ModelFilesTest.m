classdef ModelFilesTest < matlab.unittest.TestCase
    % Every live model and subsystem file loads and is saved in the team
    % release or older, so nobody pushes a model the others can't open.

    properties (Constant)
        Release = 'R2026a'
    end

    properties (TestParameter)
        file = liveModelFiles()
    end

    methods (Test)
        function notNewerThanTeamRelease(tc, file)
            info = Simulink.MDLInfo(hopperRoot(file));
            tc.verifyLessThanOrEqual(releaseOrder(info.ReleaseName), releaseOrder(tc.Release), ...
                sprintf('%s was saved in %s; export it to %s (Simulink.exportToVersion).', ...
                file, info.ReleaseName, tc.Release));
        end

        function loads(tc, file)
            [~, name] = fileparts(file);
            wasLoaded = bdIsLoaded(name);
            load_system(hopperRoot(file));
            if ~wasLoaded, tc.addTeardown(@close_system, name, 0); end
            tc.verifyTrue(bdIsLoaded(name));
        end
    end
end

function n = releaseOrder(rel)
% 'R2025b' -> 4051
n = 2 * str2double(rel(2:5)) + (lower(rel(6)) == 'b');
end

function files = liveModelFiles()
root = hopperRoot();
f = [dir(fullfile(root, 'models', '**', '*.slx')); dir(fullfile(root, 'subsystems', '**', '*.slx'))];
files = struct();
for k = 1:numel(f)
    rel = erase(fullfile(f(k).folder, f(k).name), [root filesep]);
    files.(matlab.lang.makeValidName(erase(rel, '.slx'))) = rel;
end
end

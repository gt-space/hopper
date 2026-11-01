classdef MonteCarloTest < matlab.unittest.TestCase
    % Seeded 10-run Monte Carlo (plus nominal) of the master model. Every
    % run must finish without errors or divergence and meet the CI limits
    % in mcLimits.m. Results go to outputs/montecarlo/ for CI artifacts.

    properties (Constant)
        N = 10
    end

    properties
        R
    end

    methods (TestClassSetup)
        function runAll(tc)
            outDir = hopperRoot('outputs', 'montecarlo');
            if ~isfolder(outDir), mkdir(outDir); end
            useParallel = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
            tc.R = runMonteCarlo(tc.N, UseParallel=useParallel, Limits=mcLimits(), ...
                OutputFile=fullfile(outDir, 'mc_ci.mat'));
            writetable(tc.R.summary, fullfile(outDir, 'mc_summary.csv'));
        end
    end

    methods (Test)
        function nominalPasses(tc)
            n = tc.R.nominal;
            tc.verifyTrue(n.status.success, 'nominal run errored or diverged');
            tc.verifyTrue(n.status.pass, ['nominal failed: ' strjoin(n.status.failed, ', ')]);
        end

        function everyRunCompletes(tc)
            s = tc.R.summary;
            tc.verifyTrue(all(s.success), sprintf('runs that errored or diverged: %s', mat2str(s.run(~s.success)')));
        end

        function everyRunMeetsLimits(tc)
            s = tc.R.summary;
            bad = s(~s.pass, {'run', 'failed'});
            tc.verifyEmpty(bad, sprintf('runs outside limits:\n%s', formattedDisplayText(bad)));
        end
    end
end

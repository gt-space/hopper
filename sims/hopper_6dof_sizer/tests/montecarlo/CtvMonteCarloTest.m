classdef CtvMonteCarloTest < matlab.unittest.TestCase
    % Seeded 5-run Monte Carlo (plus nominal) of the CTV on the variant
    % harness (VEHICLE=ctv, inputs/mc_params_ctv.json). The nominal flight
    % must meet the CTV's mission requirements, and every dispersed run
    % must get through setup and simulation: a run may diverge (the CTV's
    % controller isn't robust to its dispersions yet) but not error.
    % Dispersed runs are reported in outputs/montecarlo/mc_ctv_summary.csv,
    % not required to meet requirements: the CTV's geometry and tuning are
    % still provisional (docs/VARIANTS.md, The CTV).

    properties (Constant)
        N = 5
    end

    properties
        R
    end

    methods (TestClassSetup)
        function runAll(tc)
            outDir = hopperRoot('outputs', 'montecarlo');
            if ~isfolder(outDir), mkdir(outDir); end
            useParallel = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
            tc.R = runMonteCarlo(tc.N, Variants=variantConfig('VEHICLE', 'ctv'), UseParallel=useParallel, ...
                OutputFile=fullfile(outDir, 'mc_ctv_ci.mat'));
            writetable(tc.R.summary, fullfile(outDir, 'mc_ctv_summary.csv'));
        end
    end

    methods (Test)
        function nominalMeetsRequirements(tc)
            n = tc.R.nominal;
            tc.verifyTrue(n.status.success, 'nominal CTV run errored or diverged');
            tc.verifyTrue(n.status.meets_requirements, ...
                ['nominal CTV missed: ' strjoin(string(n.status.missed_requirements), ', ')]);
        end

        function noRunErrors(tc)
            % divergence is a flight outcome; anything else is a setup or code error
            s = tc.R.summary;
            bad = ~s.success & ~startsWith(s.failed, "Simulation diverged");
            tc.verifyFalse(any(bad), sprintf('CTV runs that errored:\n%s', ...
                formattedDisplayText(s(bad, {'run', 'failed'}))));
        end
    end
end

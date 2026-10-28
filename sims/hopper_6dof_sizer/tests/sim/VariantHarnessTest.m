classdef VariantHarnessTest < matlab.unittest.TestCase
    % models/hopper_6dof.slx: with every axis at its default it flies v2's
    % nominal flight (same baseline as NominalFlightTest), and every other
    % choice of every axis compiles and runs.

    properties (Constant)
        Model = 'hopper_6dof'
        SmokeTime = '2'   % s of flight per choice: compile, initialize, lift off
    end

    properties (TestParameter)
        choice = VariantHarnessTest.otherChoices()
    end

    properties
        W
    end

    methods (TestClassSetup)
        function setup(tc)
            evalin('base', 'clearvars');
            tc.W = hopperWorkspace();
        end
    end

    methods (Test)
        function defaultsMatchV2Baseline(tc)
            out = sim(buildSimInput(tc.Model, [], tc.W, variantConfig()));
            B = load(NominalFlightTest.BaselineFile);
            K = NominalFlightTest.keyNumbers(out);
            for f = fieldnames(B.K)'
                tc.verifyEqual(K.(f{1}), B.K.(f{1}), 'AbsTol', NominalFlightTest.Tol.(f{1}), f{1});
            end
        end

        function constantMassPropsHoldStill(tc)
            % MASS_PROPS=constant: no propellant flow, so the propellant
            % masses and the axial CG don't move. The CG still shifts
            % sideways by under a millimetre as the engine gimbals (its mass
            % moves with the TVC), which is physics, not propellant.
            in = buildSimInput(tc.Model, [], tc.W, variantConfig('MASS_PROPS', 'constant'));
            out = sim(in.setModelParameter('StopTime', '5'));
            range = @(s) VariantHarnessTest.range(out, s);
            tc.verifyEqual(range('ox_mass'), 0, 'ox_mass');
            tc.verifyEqual(range('fuel_mass'), 0, 'fuel_mass');
            cg = range('cg');
            tc.verifyEqual(cg(1), 0, 'axial CG');
            tc.verifyLessThan(cg(2:3), 2e-3, 'lateral CG shift from the engine gimbal');
        end

        function choiceRuns(tc, choice)
            V = variantConfig(choice{:});
            W = tc.W;
            if V.VARIANT_VEHICLE ~= 1 || V.VARIANT_CONTROLLER == 2
                W = [];   % another parameter set, or the LQI gains: build it
            end
            in = buildSimInput(tc.Model, [], W, V);
            in = in.setModelParameter('StopTime', tc.SmokeTime);
            out = sim(in);
            tc.verifyEmpty(out.ErrorMessage);
            if V.VARIANT_CONTROLLER ~= 2   % LQI is experimental: it may not hold the vehicle up
                tc.verifyEqual(out.tout(end), str2double(tc.SmokeTime), 'AbsTol', 1e-9, ...
                    sprintf('%s stopped early', strjoin(string(choice), '=')));
            end
        end
    end

    methods (Static)
        function r = range(out, signal)
            % max - min over time of each component of a logged signal
            x = squeeze(out.get(signal).Data);
            n = numel(out.tout);
            if size(x, 1) ~= n, x = reshape(permute(x, [ndims(x), 1:ndims(x)-1]), n, []); end
            r = max(x, [], 1) - min(x, [], 1);
        end

        function c = otherChoices()
            % {axis, choice} for every non-default choice, named AXIS_choice
            c = struct();
            A = variantAxes();
            for f = fieldnames(A)'
                for k = 2:numel(A.(f{1}))
                    name = char(A.(f{1})(k));
                    c.([f{1} '_' name]) = {f{1}, name};
                end
            end
            % The fluid network takes ~1 s to build thrust; without a
            % ground the vehicle falls through the pad first
            c.PROPULSION_fluid_network = [c.PROPULSION_fluid_network, {'GROUND', 'strut'}];
        end
    end
end

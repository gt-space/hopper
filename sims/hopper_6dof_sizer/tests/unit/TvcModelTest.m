classdef TvcModelTest < matlab.unittest.TestCase
    % Electromechanical TVC model (src/tvc) and its TVC_MODEL variant setup.
    % VariantHarnessTest flies the ag9113 / ag9313 choices in the harness.

    methods (Test)
        function variantAxis(tc)
            A = variantAxes();
            tc.verifyEqual(A.TVC_MODEL, ["second_order", "ag9113", "ag9313"]);
            tc.verifyEqual(variantConfig().VARIANT_TVC_MODEL, 1);
            tc.verifyEqual(variantConfig('TVC_MODEL', 'ag9313').VARIANT_TVC_MODEL, 3);
        end

        function recommendedLayoutsFitTheActuator(tc)
            % actuator at mid-stroke at 0 deg: neutral length = datasheet M
            for n = ["AG9113", "AG9313"]
                P = tvcActuatorParams(char(n), Layout='recommended');
                tc.verifyEqual(P.L0, P.M, 'AbsTol', 1e-3, n);
            end
            tc.verifyEqual(tvcActuatorParams('AG9113', Layout='recommended').r_eff, 0.1225, 'AbsTol', 1e-3);
            tc.verifyEqual(tvcActuatorParams('AG9313', Layout='recommended').r_eff, 0.0575, 'AbsTol', 1e-3);
        end

        function derivedMotorConstants(tc)
            % no-load 7/6 of rated speed, rated speed at rated current
            P = tvcActuatorParams('AG9313');
            tc.verifyEqual(P.Ke, 24 / 0.035, 'RelTol', 1e-12);
            tc.verifyEqual(P.KF, 3000 / 7.7, 'RelTol', 1e-12);
            tc.verifyEqual((P.Vbus - P.R0 * P.Imax) / P.Ke, 0.030, 'RelTol', 1e-12);
        end

        function smallStepMatchesPdr(tc)
            % PDR slide 171: AG 9113, 0.5 deg step rises in ~47 ms (10-90 %)
            S = tvcAxisSim(tvcActuatorParams('AG9113'), [], 0.6, @(t) deg2rad(0.5) * (t >= 0));
            r = TvcModelTest.riseTime(S, deg2rad(0.5));
            tc.verifyGreaterThan(r, 0.040);
            tc.verifyLessThan(r, 0.060);
            tc.verifyEqual(S.delta(end), deg2rad(0.5), 'AbsTol', deg2rad(0.02));
        end

        function largeStepIsSlewLimited(tc)
            % 8 deg: back-EMF slew; recommended AG 9113 meets 25 deg/s,
            % the sizer's AG 9313 layout does not
            A = deg2rad(8);
            S1 = tvcAxisSim(tvcActuatorParams('AG9113', Layout='recommended'), [], 0.8, @(t) A * (t >= 0));
            S2 = tvcAxisSim(tvcActuatorParams('AG9313'), [], 0.8, @(t) A * (t >= 0));
            tc.verifyGreaterThan(0.8 * 8 / TvcModelTest.riseTime(S1, A), 25);
            tc.verifyLessThan(0.8 * 8 / TvcModelTest.riseTime(S2, A), 25);
        end

        function workspaceForTheVariant(tc)
            W.IN.tvc.actuator_latency = 0.010;
            W.IN.off_axis_x = 0;  W.IN.off_axis_y = 0.02;
            W = tvcWorkspace(W, variantConfig('TVC_MODEL', 'ag9113'));
            T = W.TVCEM;
            % the Electromech choice's input Unit Delay is one of the 10 periods
            tc.verifyEqual(T.lat_steps, 9);
            tc.verifyEqual(T.e, [0.01; 0.02]);          % 1 cm floor, then the offset
            tc.verifyEqual(T.max_angle, deg2rad(8), 'AbsTol', 1e-12);
            % a TVCEM already in W is kept; Layout / Offset replace it
            W2 = tvcWorkspace(W, variantConfig('TVC_MODEL', 'ag9313'));
            tc.verifyEqual(W2.TVCEM, T);
            W3 = tvcWorkspace(W, variantConfig('TVC_MODEL', 'ag9313'), Layout='sizer', Offset=0.025);
            tc.verifyEqual(W3.TVCEM.KF, 3000 / 7.7, 'RelTol', 1e-12);
            tc.verifyEqual(W3.TVCEM.e, [0.025; 0.025]);
        end

        function updateHoldsStillAtZero(tc)
            % zero command, zero thrust: nothing moves
            T = tvcEmParams(tvcActuatorParams('AG9113'), Offset=[0; 0]);
            S = tvcEmInit();
            for k = 1:50
                [S, d] = tvcEmUpdate(S, [0; 0], 0, T);
            end
            tc.verifyEqual(d, [0; 0], 'AbsTol', 1e-12);
        end
    end

    methods (Static)
        function r = riseTime(S, A)
            d = S.delta / A;
            r = S.t(find(d >= 0.9, 1)) - S.t(find(d >= 0.1, 1));
        end
    end
end

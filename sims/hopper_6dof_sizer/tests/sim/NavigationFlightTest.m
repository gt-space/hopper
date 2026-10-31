classdef NavigationFlightTest < matlab.unittest.TestCase
    % Full hopper flights with the nav team's GMEKF (NAVIGATION axis):
    % in shadow the flight is untouched and the estimate tracks truth; on
    % the estimate the vehicle still flies the mission and lands.

    properties (Constant)
        Model = 'hopper_6dof'
        % estimate error bounds over the nominal flight (measured
        % 2026-10-29: 0.04 m, 0.13 m/s, 0.05 rad/s, 0.014), with margin
        MaxPosErr  = 0.15    % m
        MaxVelErr  = 0.4     % m/s
        MaxRateErr = 0.15    % rad/s
        MaxQuatErr = 0.04    % |q_est - q_true|
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
        function shadowTracksTruth(tc)
            out = sim(buildSimInput(tc.Model, [], tc.W, variantConfig('NAVIGATION', 'gmekf_shadow')));
            % the controller still flies on truth: v2's flight
            B = load(NominalFlightTest.BaselineFile);
            K = NominalFlightTest.keyNumbers(out);
            for f = fieldnames(B.K)'
                tc.verifyEqual(K.(f{1}), B.K.(f{1}), 'AbsTol', NominalFlightTest.Tol.(f{1}), f{1});
            end
            tc.verifyEstimate(out);
        end

        function flyingOnTheEstimateLands(tc)
            out = sim(buildSimInput(tc.Model, [], tc.W, variantConfig('NAVIGATION', 'gmekf')));
            K = NominalFlightTest.keyNumbers(out);
            IN = tc.W.IN;
            tc.verifyGreaterThan(K.max_altitude, IN.mission.target_altitude, 'apex');
            tc.verifyLessThan(abs(K.final_altitude), 0.1, 'landed');
            tc.verifyLessThan(K.landing_offset, IN.mission.max_landing_offset, 'landing offset');
            tc.verifyLessThan(K.max_tilt_deg, IN.mission.max_tilt_deg, 'tilt');
            tc.verifyEstimate(out);
        end

        function fdirEstimatorChecksCatchGpsBias(tc)
            % FDIR on the shadow filter: quiet until a 5 m GPS bias at 10 s,
            % then the NIS CUSUM (row 35, R1) latches; the vehicle lands
            W = tc.W;
            W.faults = faultScenario({'gps', 'bias', 10, inf, 5});
            V = variantConfig('FDIR', 'on', 'NAVIGATION', 'gmekf_shadow', 'FAULT_INJECTION', 'on');
            out = sim(buildSimInput(tc.Model, [], W, V));
            st = NavigationFlightTest.logged(out, 'fdir_status');
            t = out.get('fdir_status').Time;
            tc.verifyFalse(any(any(st(t < 10, 33:37) == 2)), 'estimator fault before the GPS bias');
            tc.verifyTrue(any(st(:, 35) == 2), 'NIS CUSUM (row 35) never latched');
            M = fdir_modes();
            mode = NavigationFlightTest.logged(out, 'fdir_mode');
            tc.verifyEqual(mode(end), double(M.LANDED), 'final mode');
        end

        function fdirOnTheEstimateFliesQuietly(tc)
            % FDIR monitors on the estimate (NAVIGATION gmekf): no fault, PDR
            % mode sequence to LANDED
            V = variantConfig('FDIR', 'on', 'NAVIGATION', 'gmekf');
            out = sim(buildSimInput(tc.Model, [], tc.W, V));
            tc.verifyEqual(max(NavigationFlightTest.logged(out, 'fdir_level')), -1, 'response level');
            M = fdir_modes();
            mode = NavigationFlightTest.logged(out, 'fdir_mode');
            tc.verifyEqual(mode(end), double(M.LANDED), 'final mode');
        end

        function ctvFilterRuns(tc)
            % The CTV's GMEKF compiles and runs in shadow. Only a smoke run:
            % the CTV plant spins up in roll from ~0.8 s (its roll loop, not
            % the filter), and no estimator can follow that yet
            V = variantConfig('VEHICLE', 'ctv', 'NAVIGATION', 'gmekf_shadow');
            in = buildSimInput(tc.Model, [], [], V);
            out = sim(in.setModelParameter('StopTime', '0.6'));
            tc.verifyEmpty(out.ErrorMessage);
            E = NavigationFlightTest.logged(out, 'nav_error');
            tc.verifyLessThan(max(vecnorm(E(:, 1:3), 2, 2)), tc.MaxPosErr, 'position error (m)');
        end
    end

    methods
        function verifyEstimate(tc, out)
            E = NavigationFlightTest.logged(out, 'nav_error');
            tc.verifyLessThan(max(vecnorm(E(:, 1:3), 2, 2)), tc.MaxPosErr, 'position error (m)');
            tc.verifyLessThan(max(vecnorm(E(:, 4:6), 2, 2)), tc.MaxVelErr, 'velocity error (m/s)');
            tc.verifyLessThan(max(vecnorm(E(:, 7:9), 2, 2)), tc.MaxRateErr, 'rate error (rad/s)');
            tc.verifyLessThan(max(vecnorm(E(:, 10:13), 2, 2)), tc.MaxQuatErr, 'attitude error');
            Pd = NavigationFlightTest.logged(out, 'nav_P_diag');
            tc.verifyTrue(all(isfinite(Pd(:))) && all(Pd(:) > 0), 'covariance stays positive');
        end
    end

    methods (Static)
        function x = logged(out, name)
            s = out.get(name);
            x = double(squeeze(s.Data));
            if size(x, 1) ~= numel(s.Time), x = x'; end
        end
    end
end

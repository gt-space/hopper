classdef NominalFlightTest < matlab.unittest.TestCase
    % Nominal hopper flight through buildSimInput, against stored key
    % numbers. Regenerate the baseline after an intended physics change with
    %   NominalFlightTest.updateBaseline()

    properties (Constant)
        Model = 'hopper_6dof_NED_v2'
        BaselineFile = hopperRoot('tests', 'baselines', 'nominal_v2.mat')
        % Absolute tolerance per key number. The same flight on Linux (CI)
        % and Windows differs by about a millimetre at touchdown (different
        % math libraries), so bit-exactness can't hold across platforms.
        % These are still far below any real change to the flight.
        Tol = struct('flight_time', 0.01, ...     % s
                     'max_altitude', 0.01, ...    % m
                     'final_altitude', 0.01, ...  % m
                     'landing_offset', 0.02, ...  % m
                     'final_velocity', 0.02, ...  % m/s
                     'max_tilt_deg', 0.05, ...    % deg
                     'final_ox_mass', 0.01, ...   % kg
                     'final_fu_mass', 0.01)       % kg
    end

    properties
        out
        W
    end

    methods (TestClassSetup)
        function fly(tc)
            evalin('base', 'clearvars');
            [in, tc.W] = buildSimInput(tc.Model);
            tc.out = sim(in);
        end
    end

    methods (Test)
        function runsWithoutBaseWorkspace(tc)
            tc.verifyEmpty(tc.out.ErrorMessage);
            tc.verifyEmpty(evalin('base', 'who'), 'buildSimInput must not leave variables in base');
        end

        function matchesBaseline(tc)
            B = load(tc.BaselineFile);
            K = NominalFlightTest.keyNumbers(tc.out);
            for f = fieldnames(B.K)'
                tc.verifyEqual(K.(f{1}), B.K.(f{1}), 'AbsTol', tc.Tol.(f{1}), ...
                    sprintf('%s: got %s, baseline %s (tolerance %g)', f{1}, ...
                    mat2str(K.(f{1}), 6), mat2str(B.K.(f{1}), 6), tc.Tol.(f{1})));
            end
        end

        function reachesTargetAndLands(tc)
            K = NominalFlightTest.keyNumbers(tc.out);
            tc.verifyGreaterThan(K.max_altitude, 51);
            tc.verifyLessThan(abs(K.final_altitude), 0.1);
        end

        function sloshOpenLoopMatchesTruth(tc)
            % SLOSH_MODEL.md section 7.3: Slosh subsystem vs Newton-Euler truth
            for v = {'OUT', 'MoI_init', 'cg_init', 'TANKS'}
                assignin('base', v{1}, tc.W.(v{1}));
            end
            tc.addTeardown(@() evalin('base', 'clearvars'));
            cases = [1 0; 0.1 0; 0 0.25; 0 1];   % LOx 100/10 %, RP-1 25/100 %
            for k = 1:size(cases, 1)
                R = sloshOpenLoopTest([tc.Model '.slx'], cases(k, 1), cases(k, 2), 'specific', 'asbuilt');
                label = sprintf('ox %.2f / fu %.2f', cases(k, 1), cases(k, 2));
                tc.verifyLessThan(R.err_th, 0.01, ['pitch error, ' label]);
                tc.verifyLessThan(R.err_y, 0.01, ['bob displacement error, ' label]);
                tc.verifyGreaterThan(R.F0sign, 0.9, ['initial force sign, ' label]);
                tc.verifyLessThan(R.F0sign, 1.2, ['initial force sign, ' label]);
            end
        end
    end

    methods (Static)
        function K = keyNumbers(out)
            z = squeeze(out.get('z').Data);
            x = squeeze(out.get('x').Data); y = squeeze(out.get('y').Data);
            v = squeeze(out.get('velocity').Data);
            p = squeeze(out.get('pitch').Data);
            K.flight_time    = out.tout(end);
            K.max_altitude   = -min(z);
            K.final_altitude = -z(end);
            K.landing_offset = hypot(x(end), y(end));
            K.final_velocity = v(end, :);
            K.max_tilt_deg   = max(abs(90 - p));
            K.final_ox_mass  = out.get('ox_mass').Data(end);
            K.final_fu_mass  = out.get('fuel_mass').Data(end);
        end

        function updateBaseline()
            evalin('base', 'clearvars');
            out = sim(buildSimInput(NominalFlightTest.Model)); %#ok<PROP>
            K = NominalFlightTest.keyNumbers(out);
            K %#ok<NOPRT>
            save(NominalFlightTest.BaselineFile, 'K');
        end
    end
end

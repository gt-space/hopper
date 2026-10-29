classdef ModeManagerTest < matlab.unittest.TestCase
    % fdir_mode_manager / fdir_touchdown against the GNC PDR state machine
    % (slides 10-11) and contingency rules (slides 79, 84-86). No Simulink.

    properties
        M
        P
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.M = fdir_modes();
            IN.mission.trajectory = struct('apex', 52, 't_hover', 2, 'z_td', 0.4);
            tc.P = fdir_mode_params(IN);
        end
    end

    methods (Test)
        function nominalFlightVisitsEveryMode(tc)
            M = tc.M;
            s = fdir_mode_state();
            seen = fly(tc, s, nominalProfile());
            expected = [M.SAFE M.INITIALIZE M.READY M.ARMED M.ENGINE_START M.PAD_CHECKOUT ...
                M.ASCENT M.HOVER M.DESCENT M.TERMINAL_DESCENT M.TOUCHDOWN_CONFIRM M.LANDED];
            tc.verifyEqual(unique(seen, 'stable'), expected);
        end

        function modeChangeTakesEffectNextFrame(tc)
            s = fdir_mode_state();
            [out1, s] = fdir_mode_manager(frame(0, 0, 0), tc.P, s);
            out2 = fdir_mode_manager(frame(0.01, 0, 0), tc.P, s);
            tc.verifyEqual(out1.mode, tc.M.SAFE);
            tc.verifyEqual(out2.mode, tc.M.INITIALIZE);
        end

        function r2InFlightAbortsToLand(tc)
            s = inMode(tc, tc.M.ASCENT);
            u = frame(5, 20, 3); u.response = tc.M.R2;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.ABORT_LAND);
        end

        function r0AndR1NeverChangeTheMode(tc)
            s = inMode(tc, tc.M.HOVER);
            for r = [tc.M.R0 tc.M.R1]
                u = frame(12.5, 52, 0); u.response = r;   % inside the 2 s hover
                s = step2(tc, u, s);
                tc.verifyEqual(s.mode, tc.M.HOVER);
            end
        end

        function r4AboveSafeHeightFliesAbortFirst(tc)
            % automatic R4 in flight waits for h_cut (PDR slide 79)
            s = inMode(tc, tc.M.ASCENT);
            u = frame(5, 20, 3); u.response = tc.M.R4;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.ABORT_LAND);
            tc.verifyTrue(s.r4_pending);
            u = frame(30, tc.P.h_cut - 0.1, -0.5); u.response = tc.M.R4;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.SAFED);
        end

        function r4WhenUncontrollableSafesAtAnyHeight(tc)
            s = inMode(tc, tc.M.ASCENT);
            u = frame(5, 20, 3); u.response = tc.M.R4; u.uncontrollable = true;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.SAFED);
        end

        function terminateIsAlwaysAccepted(tc)
            s = inMode(tc, tc.M.HOVER);
            u = frame(14, 52, 0); u.cmd.terminate = true;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.SAFED);
        end

        function r3IsNavDegradedAbort(tc)
            s = inMode(tc, tc.M.DESCENT);
            u = frame(20, 30, -3); u.response = tc.M.R3;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.NAV_ABORT);
        end

        function abortOnThePadIsAShutdown(tc)
            s = inMode(tc, tc.M.ENGINE_START);
            u = frame(0.1, 0, 0); u.response = tc.M.R2;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.SAFED);
        end

        function noArmHoldsInReady(tc)
            s = inMode(tc, tc.M.READY);
            u = frame(0.05, 0, 0); u.no_arm = true;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.READY);
        end

        function abortExitsOnlyOnTouchdown(tc)
            s = inMode(tc, tc.M.ABORT_LAND);
            s = step2(tc, frame(25, 0.5, -0.3), s);
            tc.verifyEqual(s.mode, tc.M.ABORT_LAND);
            u = frame(26, 0, 0); u.touchdown = true;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.LANDED);
        end

        function groundResetFromLanded(tc)
            s = inMode(tc, tc.M.LANDED);
            u = frame(40, 0, 0); u.cmd.reset = true;
            s = step2(tc, u, s);
            tc.verifyEqual(s.mode, tc.M.SAFE);
        end

        function touchdownNeedsTwoSignalsAndPersistence(tc)
            M = tc.M; P = tc.P;
            s = fdir_mode_state();
            one = struct('lidar_median', 0.01, 'accel_impulse', 0, 'vz', -0.5, 'thrust_cut', false);
            for k = 1:10
                [c, n, s] = fdir_touchdown(one, P, M.TERMINAL_DESCENT, s);
            end
            tc.verifyEqual(n, uint8(1));
            tc.verifyFalse(c, 'one signal alone must not confirm');
            two = one; two.thrust_cut = true;
            confirmedAt = 0;
            for k = 1:10
                [c, ~, s] = fdir_touchdown(two, P, M.TERMINAL_DESCENT, s);
                if c && ~confirmedAt, confirmedAt = k; end
            end
            tc.verifyEqual(confirmedAt, P.td_n, 'confirms on the N-th frame with two signals');
        end

        function touchdownDisarmedAtAltitude(tc)
            s = fdir_mode_state();
            all4 = struct('lidar_median', 0, 'accel_impulse', 10, 'vz', 0, 'thrust_cut', true);
            for k = 1:10
                [c, ~, s] = fdir_touchdown(all4, tc.P, tc.M.HOVER, s);
            end
            tc.verifyFalse(c);
        end
    end

    methods
        function seen = fly(tc, s, prof)
            seen = zeros(1, numel(prof), 'uint8');
            for k = 1:numel(prof)
                [out, s] = fdir_mode_manager(prof(k), tc.P, s);
                seen(k) = out.mode;
            end
        end

        function s = inMode(~, mode)
            s = fdir_mode_state();
            s.mode = mode;
            s.hover_start = 12;
        end

        function s = step2(tc, u, s)
            [~, s] = fdir_mode_manager(u, tc.P, s);
        end
    end
end

function u = frame(t, alt, vz)
u.t = t; u.alt = alt; u.vz = vz;
u.thrust_ratio = 1.1; u.response = uint8(0); u.uncontrollable = false;
u.touchdown = false; u.thrust_cut = false; u.no_arm = false;
u.cmd = struct('init', false, 'arm', false, 'start', false, 'launch', false, ...
    'abort', false, 'terminate', false, 'reset', false);
end

function prof = nominalProfile()
% 52 m hop: pad, climb, 2 s hover, descent, terminal, cut, settled
k = 0;
add = @(t, alt, vz) frame(t, alt, vz);
prof = add(0, 0, 0);
for t = 0.01:0.01:0.1, k = k + 1; prof(end+1) = add(t, 0, 0); end %#ok<AGROW>
prof(end+1) = add(0.2, 0.05, 0.5);
for t = 1:12, prof(end+1) = add(t, 52 * t / 12, 2); end
prof(end+1) = add(12.1, 52, 0.1);
for t = 12.2:0.5:14.6, prof(end+1) = add(t, 52, 0); end
for t = 15:1:27, prof(end+1) = add(t, max(0.3, 52 - 4 * (t - 14)), -3); end
u = add(28, 0.04, -0.3); u.thrust_cut = true; prof(end+1) = u;
for t = 28.1:0.1:28.6
    u = add(t, 0, 0); u.thrust_cut = true; u.touchdown = t > 28.3; prof(end+1) = u; %#ok<AGROW>
end
end

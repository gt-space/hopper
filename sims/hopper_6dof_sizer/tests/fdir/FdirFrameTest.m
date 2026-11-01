classdef FdirFrameTest < matlab.unittest.TestCase
    % fdir_frame end to end without Simulink: a steady hover with no fault
    % stays quiet; injected faults latch the right catalog rows and drive
    % the response ladder (FDIR CoDR 15-28, GNC PDR 78-85).

    properties
        P
        M
        K
    end

    methods (TestMethodSetup)
        function setup(tc)
            IN = mission_inputs();
            S = sim_sensor_params(IN);
            tc.P.S = S;
            tc.P.mon = gnc_monitor_params(IN, S);
            tc.P.modes = fdir_mode_params(IN);
            tc.P.modes.contact_alt = S.lidar.mount_x / cosd(S.lidar.cant) + 0.05;
            tc.P.cat = fdir_load_catalog();
            tc.P.h_cut = tc.P.modes.h_cut;
            tc.P.abort = struct('v_land', 0.3, 'v_max', 3, 'tau', 4, 'a_ramp', 1, 'dt', S.frame);
            tc.P.faults = faultScenario();
            tc.P.inject = true;
            tc.P.ref0_err = 0;
            tc.P.thrust_max = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(2);
            tc.P.thrust_min = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(1);
            tc.P.g = 9.80665;
            tc.P.cutoff_alt = IN.mission.engine_cutoff_alt;
            tc.P.nav_active = false;                    % estimator checks: NavigationFlightTest
            tc.M = fdir_modes();
            tc.K = fdir_catalog_cols();
        end
    end

    methods (Test)
        function quietHoverLatchesNothing(tc)
            [ys, ~] = hover(tc, 3.5);                  % hover lasts 2 s from t = 2
            y = ys(end);
            tc.verifyEqual(y.mode, double(tc.M.HOVER), 'climb then hover');
            latched = find(y.status == 2);
            tc.verifyEmpty(latched, sprintf('false alarms: %s', mat2str(latched')));
            tc.verifyEqual(y.level, -1);
        end

        function lidarBeamLossMasksThatBeam(tc)
            tc.P.modes.hover_alt = 5;                   % inside the lidar's rated range
            tc.P.faults = faultScenario({'lidar', 'loss', 2.5, inf, 0, 2});
            ys = hover(tc, 3.5);
            y = ys(end);
            tc.verifyEqual(y.status(9), 2, 'lidar beam 2 latched');
            tc.verifyEqual(y.level, 0, 'R0 only');
            tc.verifyTrue(bitand(y.mask, 2^(5 - 1)) > 0, 'beam 2 masked');
            tc.verifyEqual(y.mode, double(tc.M.HOVER), 'R0 never changes the mode');
        end

        function gyroStuckIsR4(tc)
            tc.P.faults = faultScenario({'gyro', 'stuck', 4});
            ys = hover(tc, 6);
            y = ys(end);
            tc.verifyEqual(y.status(5), 2, 'gyro stuck latched');
            tc.verifyEqual(y.level, 4);
            tc.verifyTrue(y.abort || y.safed, 'R4 above h_cut flies the abort until permitted');
        end

        function imuLossLatchesWithinTwoFrames(tc)
            tc.P.faults = faultScenario({'accel', 'loss', 4});
            ys = hover(tc, 4.05);
            tc.verifyEqual(ys(end).status(1), 2);
        end
    end

    methods
        function [ys, st] = hover(tc, t_end)
            % climb from the pad to the hover altitude, then hold
            % the shared Sensors block (sensor_frame), run here at the FDIR frame
            st = fdir_frame_state(size(tc.P.cat, 1), size(tc.P.faults, 1), tc.P.S.frame);
            sp = struct('S', tc.P.S, 'faults', tc.P.faults, 'inject', tc.P.inject);
            sst = sensor_frame_state(size(tc.P.faults, 1), tc.P.S.frame);
            L = sensor_layout();
            h = tc.P.modes.hover_alt;
            ts = 0:tc.P.S.frame:t_end;
            ys = cell(1, 0);
            rng(1);
            noise = randn(L.N_MEAS, numel(ts));
            for k = 1:numel(ts)
                t = ts(k);
                alt = min(h, h * t / 2);
                vz_up = (t < 2) * h / 2;
                in.t = t;
                in.x = [0; 0; -alt; 0; 0; -vz_up; 0; 0; 0; sqrt(2)/2; 0; sqrt(2)/2; 0];
                in.dcm_be = [0 0 -1; 0 1 0; 1 0 0];     % body x up
                in.f_b = [9.80665; 0; 0];
                [in.meas, in.valid, in.stamp, ~, sst] = sensor_frame(in.x, in.f_b, in.dcm_be, [1.14; 0; 0], ...
                    noise(:, k), t, sp, sst);
                in.thrust = 1300; in.gimbal_fb = [0; 0]; in.rcs_torque = 0;
                in.mass = 130; in.u_nom = [1300; 0; 0; 0]; in.u_cmd = in.u_nom;
                [y, st] = fdir_frame(in, tc.P, st);
                ys{end + 1} = y; %#ok<AGROW>
            end
            ys = [ys{:}];
        end
    end
end

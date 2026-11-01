classdef NavAdapterTest < matlab.unittest.TestCase
    % The adapter between the sim and the nav team's GMEKF (src/navigation):
    % frames and states convert exactly both ways, and with noise off each
    % sensor's measurement from our sensor model is what the filter's own
    % measurement model predicts at the true state.

    properties
        F      % nav_frames for the hopper on the pad
        cg     % a CG away from the pad value, as mid-flight
    end

    methods (TestMethodSetup)
        function frames(tc)
            dcm_pad = angle2dcm(0, pi/2, 0);    % the 6DOF block's eul_0 = [0 pi/2 0]
            tc.F = nav_frames(dcm_pad, 1.1429);
            tc.cg = [1.08; 0.002; -0.001];
        end
    end

    methods (Test)
        function padIsNavIdentity(tc)
            % on the pad the nav body axes are N, E, D (their hover tests)
            tc.verifyEqual(tc.F.R * angle2dcm(0, pi/2, 0), eye(3), 'AbsTol', 1e-15);
            tc.verifyEqual(tc.F.R * [1; 0; 0], [0; 0; -1], 'AbsTol', 1e-15, 'sim +x (up) is nav -z');
        end

        function quaternionConventions(tc)
            % the sim's quat2dcm and the nav's quaternion_to_rotation' are the
            % same A(q), and nav_quat_from_dcm inverts it
            rng(3);
            for k = 1:20
                q = randn(4, 1); q = q / norm(q); q = q * sign(q(1));
                tc.verifyEqual(nav_dcm_from_quat(q), quat2dcm(q'), 'AbsTol', 1e-14);
                tc.verifyEqual(nav_quat_from_dcm(quat2dcm(q')), q, 'AbsTol', 1e-12);
            end
        end

        function stateRoundTrip(tc)
            rng(4);
            for k = 1:20
                eul = [0 pi/2 0] + deg2rad(20) * (rand(1, 3) - 0.5);
                C = angle2dcm(eul(3), eul(2), eul(1));
                x = [randn(3, 1) * 10; randn(3, 1); randn(3, 1) * 0.2; nav_quat_from_dcm(C)];
                xn = nav_state_from_sim(x, C, tc.cg, tc.F);
                [x2, C2] = nav_state_to_sim(xn, tc.cg, tc.F);
                tc.verifyEqual(x2, x, 'AbsTol', 1e-12);
                tc.verifyEqual(C2, C, 'AbsTol', 1e-12);
            end
        end

        function estimateKeepsTheFilterSign(tc)
            % the controller's quaternion follows the filter's sign (q and
            % -q are one attitude): no flip where w crosses zero
            rng(6);
            for k = 1:20
                q = randn(4, 1); q = q / norm(q);
                xn = [zeros(9, 1); q; zeros(7, 1)];
                [x1, C1] = nav_state_to_sim(xn, tc.cg, tc.F);
                xn(10:13) = -q;
                [x2, C2] = nav_state_to_sim(xn, tc.cg, tc.F);
                tc.verifyEqual(x2(10:13), -x1(10:13), 'AbsTol', 1e-15, 'sign follows the filter');
                tc.verifyEqual(C2, C1, 'AbsTol', 1e-15, 'same attitude');
                tc.verifyEqual(nav_dcm_from_quat(x1(10:13)), C1, 'AbsTol', 1e-12, 'quaternion matches the DCM');
            end
        end

        function measurementsMatchFilterModels(tc)
            % noise off: at the true state, each of the nav team's
            % measurement models predicts what our sensor model measured
            IN = mission_inputs();
            S = sim_sensor_params(IN);
            P.F = tc.F; P.S = S;
            [P.geom.lidar_dir, P.geom.lidar_mount] = lidar_geometry(S.lidar);
            P.geom.gps_antenna = S.gps.antenna;
            L = sensor_layout();
            rng(5);
            for k = 1:10
                eul = [0 pi/2 0] + deg2rad(10) * (rand(1, 3) - 0.5);
                C = angle2dcm(eul(3), eul(2), eul(1));
                x = [randn(2, 1); -5 - 10 * rand; randn(3, 1); 0.1 * randn(3, 1); nav_quat_from_dcm(C)];
                xs = struct('pos_ned', x(1:3), 'vel_ned', x(4:6), 'w_b', x(7:9), 'dcm_be', C, ...
                            'f_b', [9.8; 0; 0], 'cg', tc.cg);
                [~, tr] = sensor_models(xs, zeros(L.N_MEAS, 1), 0, S, sensor_state());
                [~, angvel, ~, gps, lidar, dgps, mag, lf, gf] = nav_measurements(tr, true(L.N_SENSORS, 1), ...
                    zeros(L.N_SENSORS, 1), tc.cg, 130, eye(3), zeros(3, 1), P, 0, 1);
                xn = nav_state_from_sim(x, C, tc.cg, tc.F);
                inn = gyro_update_model(xn, angvel, struct());
                tc.verifyEqual(inn, zeros(3, 1), 'AbsTol', 1e-12, 'gyro');
                inn = gps_update_model(xn, gps, struct('r_b', gps_lever_arm(gf)));
                tc.verifyEqual(inn, zeros(6, 1), 'AbsTol', 1e-9, 'GPS position and velocity');
                inn = diffGps_update_model(xn, dgps, gf);
                tc.verifyEqual(inn, zeros(3, 1), 'AbsTol', 1e-12, 'differential GPS');
                inn = mag_update_model(xn, mag, struct());
                tc.verifyEqual(inn, zeros(3, 1), 'AbsTol', 1e-12, 'magnetometer');
                for b = 1:4
                    lp = lf; lp.lidarIndex = b;
                    inn = lidar_update_model(xn, struct('distance', lidar.distance(b)), lp);
                    tc.verifyEqual(inn, 0, 'AbsTol', 1e-9, sprintf('lidar %d', b));
                end
            end
        end

        function fdirMasksReachTheFilter(tc)
            % a sensor FDIR masks goes to the filter as invalid; the R scale
            % inflates every covariance; the mask bits are fdir_targets'
            IN = mission_inputs();
            P.F = tc.F; P.S = sim_sensor_params(IN);
            [P.geom.lidar_dir, P.geom.lidar_mount] = lidar_geometry(P.S.lidar);
            P.geom.gps_antenna = P.S.gps.antenna;
            L = sensor_layout();
            meas = zeros(L.N_MEAS, 1); meas(L.LIDAR) = 2;
            args = {meas, true(L.N_SENSORS, 1), zeros(L.N_SENSORS, 1), tc.cg, 130, eye(3), zeros(3, 1), P};
            T = fdir_targets();
            bit = @(names) sum(arrayfun(@(n) 2^(find(T == n) - 1), string(names)));
            [a0, w0, b0, g0, l0, d0] = nav_measurements(args{:}, 0, 1);
            tc.verifyTrue(a0.validity && w0.validity && b0.validity && g0.validity && all(l0.validity) && d0.validity);
            [a, w, b, g, l, d] = nav_measurements(args{:}, bit(["imu_gyro", "lidar3", "gps"]), 4);
            tc.verifyTrue(a.validity, 'accel stays');
            tc.verifyFalse(w.validity, 'gyro masked');
            tc.verifyEqual(l.validity(:)', [true true false true], 'lidar 3 masked');
            tc.verifyFalse(g.validity || d.validity, 'GPS and differential GPS masked');
            tc.verifyTrue(b.validity, 'baro stays');
            tc.verifyEqual(b.covariance, 4 * b0.covariance, 'AbsTol', 1e-15, 'R scale');
            [a, w] = nav_measurements(args{:}, bit("imu"), 1);
            tc.verifyFalse(a.validity || w.validity, 'IMU target masks both');
        end

        function padStateIsNavZero(tc)
            % standing on the pad, CG at its rest height: nav position 0,
            % identity attitude
            C = angle2dcm(0, pi/2, 0);
            cg0 = [tc.F.heightCOM; 0; 0];
            x = [sim_pos_offset(cg0) + [0; 0; -tc.F.heightCOM]; zeros(3, 1); zeros(3, 1); nav_quat_from_dcm(C)];
            xn = nav_state_from_sim(x, C, cg0, tc.F);
            tc.verifyEqual(xn(1:3), zeros(3, 1), 'AbsTol', 1e-15);
            tc.verifyEqual(xn(10:13), [1; 0; 0; 0], 'AbsTol', 1e-15);
        end
    end
end

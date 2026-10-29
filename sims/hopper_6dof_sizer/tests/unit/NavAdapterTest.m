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
                [~, angvel, ~, gps, lidar, dgps, ~, lf, gf] = nav_measurements(tr, true(L.N_SENSORS, 1), ...
                    zeros(L.N_SENSORS, 1), tc.cg, 130, eye(3), zeros(3, 1), P);
                xn = nav_state_from_sim(x, C, tc.cg, tc.F);
                inn = gyro_update_model(xn, angvel, struct());
                tc.verifyEqual(inn, zeros(3, 1), 'AbsTol', 1e-12, 'gyro');
                inn = gps_update_model(xn, gps, struct('r_b', gps_lever_arm(gf)));
                tc.verifyEqual(inn, zeros(6, 1), 'AbsTol', 1e-9, 'GPS position and velocity');
                inn = diffGps_update_model(xn, dgps, gf);
                tc.verifyEqual(inn, zeros(3, 1), 'AbsTol', 1e-12, 'differential GPS');
                for b = 1:4
                    lp = lf; lp.lidarIndex = b;
                    inn = lidar_update_model(xn, struct('distance', lidar.distance(b)), lp);
                    tc.verifyEqual(inn, 0, 'AbsTol', 1e-9, sprintf('lidar %d', b));
                end
            end
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

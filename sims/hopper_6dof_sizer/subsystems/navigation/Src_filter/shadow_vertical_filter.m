%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shadow_vertical_filter.m
%
% The shadow vertical estimator (GNC PDR, R3 nav-degraded abort): a small
% Kalman filter on height and vertical velocity alone, independent of the
% GMEKF's position, velocity and biases. FDIR compares it with the main
% estimate (fault table row 36) and flies the nav-degraded abort on it.
%
% State xs = [z; vz]: the COM's NED z (from its pad rest point, the GMEKF's
% convention) and its rate. Propagation uses the IMU's specific force
% rotated by the main filter's attitude (the one thing it borrows: at the
% hopper's tilts a degree of attitude error moves the vertical
% acceleration by under 0.02 %); the barometer and the LIDAR beams correct
% it with the same measurement models the GMEKF uses.
%
% Inputs:
%   xs, Ps: shadow state (2x1) and covariance (2x2)
%   accel_meas: accelBus (body specific force)
%   q_main: GMEKF attitude quaternion (4x1), body to inertial
%   meas_baro: baroBus; meas_LID: lidarBus; lidar_frames: lidarFrameBus
%   g: gravity (scalar); dt: step (s); q_accel: vertical acceleration
%   noise PSD ((m/s^2)^2/Hz)
%   st: .baro_t .lid_t timestamps already used (start at -1)
%
% Outputs:
%   xs, Ps, st updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xs, Ps, st] = shadow_vertical_filter(xs, Ps, accel_meas, q_main, meas_baro, meas_LID, ...
    lidar_frames, g, dt, q_accel, st)

    % ---- Propagate ----
    C_IB = quaternion_to_rotation(q_main);
    a_z = 0;
    if accel_meas.validity
        a_I = C_IB * accel_meas.accel(:);
        a_z = a_I(3) + g;
    end
    F = [1 dt; 0 1];
    xs = F * xs + [0.5 * dt^2; dt] * a_z;
    G = [0.5 * dt^2; dt];
    Ps = F * Ps * F' + G * G' * q_accel / dt;

    % ---- Barometer: height = -z (bias-free here) ----
    if meas_baro.validity && meas_baro.timestamp > st.baro_t
        H = [-1 0];
        [xs, Ps] = scalarUpdate(xs, Ps, H, meas_baro.height - (-xs(1)), meas_baro.covariance);
        st.baro_t = meas_baro.timestamp;
    end

    % ---- LIDAR beams: distance = -(n r + z - heightCOM) / (n A(:,1)) ----
    if meas_LID.timestamp > st.lid_t
        n = C_IB(3, :);
        for k = 1:4
            if ~meas_LID.validity(k), continue; end
            switch k
                case 1, A = lidar_frames.A1; r = lidar_frames.r1;
                case 2, A = lidar_frames.A2; r = lidar_frames.r2;
                case 3, A = lidar_frames.A3; r = lidar_frames.r3;
                otherwise, A = lidar_frames.A4; r = lidar_frames.r4;
            end
            D = n * A(:, 1);
            if abs(D) < 1e-3, continue; end
            d_pred = -(n * r + xs(1) - lidar_frames.heightCOM) / D;
            if d_pred <= 0, continue; end
            H = [-1 / D 0];
            [xs, Ps] = scalarUpdate(xs, Ps, H, meas_LID.distance(k) - d_pred, meas_LID.covariance(k));
        end
        st.lid_t = meas_LID.timestamp;
    end
end

function [x, P] = scalarUpdate(x, P, H, r, R)
% one scalar measurement, 99 % gate (chi-square, 1 DOF)
S = H * P * H' + R;
if r^2 / S > 6.635, return; end
K = P * H' / S;
x = x + K * r;
I_KH = eye(2) - K * H;
P = I_KH * P * I_KH' + K * R * K';
end

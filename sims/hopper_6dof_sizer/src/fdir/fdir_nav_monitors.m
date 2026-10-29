function [raw, metric, nis, m] = fdir_nav_monitors(meas, valid, stamp, xn, Pn, cg, c, m)
% FDIR_NAV_MONITORS  Estimator health checks, catalog rows 33-37.
%
%   meas, valid, stamp  the shared Sensors block this frame (sensor_layout)
%   xn, Pn              nav estimate (20) and covariance (19 x 19)
%   cg                  Mass_Properties CG (3), for the lever arms
%   c                   .nav (nav_params: frames, sigmas, geometry), the
%                       gates and plausibility bounds (gnc_monitor_params)
%   m                   state: .stamp (sensor stamps last checked)
%
%   raw     37 x 1 logical, rows 33, 34, 35 (unused: CUSUM), 37 set
%   metric  37 x 1, row 35's CUSUM increment
%   nis     5 x 1 normalized NIS this frame (baro gyro gps lidar dgps;
%           -1 when that sensor had no new sample), for logs and tuning
%
% Each new sample is checked against the estimate with the nav team's own
% measurement models (subsystems/navigation/Src_filter) and 99.99 % gates
% (gnc_monitor_params): NIS = r' S^-1 r, S = H P H' + R. The filter gates on the prior;
% this checks the posterior, so a sensor the filter just accepted passes
% and one it rejected (or that dragged the estimate away) does not.
%   33  exactly one sensor group fails its gate (lidar: two beams at once)
%   34  two or more fail together (divergence)
%   35  CUSUM of the mean normalized NIS above 2 (poor tuning), per second
%   37  covariance or state implausible: not finite, a variance <= 0 or
%       above its bound, quaternion off unit length, speed or rate beyond
%       physical limits
% Row 36 (main vs shadow estimator) needs a second estimator.
%
%#codegen

raw = false(37, 1); metric = zeros(37, 1); nis = -ones(5, 1);
L = sensor_layout();
xn = xn(:);
[~, angvel, baro, gps, lidar, dgps, ~, lf, gf] = nav_measurements(meas, valid, stamp, cg, 1, eye(3), zeros(3, 1), c.nav);

fail = false(5, 1);
fresh = stamp(:) ~= m.stamp(:);
if fresh(L.S_BARO) && baro.validity
    [r, H] = baro_update_model(xn, baro, struct());
    [nis(1), fail(1)] = check(r, H, Pn, baro.covariance, c.gate_baro);
end
if fresh(L.S_IMU) && angvel.validity
    [r, H] = gyro_update_model(xn, angvel, struct());
    [nis(2), fail(2)] = check(r, H, Pn, angvel.covariance, c.gate_gyro);
end
if fresh(L.S_GPS) && gps.validity
    [r, H] = gps_update_model(xn, gps, struct('r_b', gf.r1));
    [nis(3), fail(3)] = check(r, H, Pn, blkdiag(gps.pcovariance, gps.vcovariance), c.gate_gps);
end
if any(fresh(L.S_LIDAR))
    worst = -1; nbad = 0;
    for b = 1:4
        if lidar.validity(b)
            lp = lf; lp.lidarIndex = b;
            [r, H] = lidar_update_model(xn, struct('distance', lidar.distance(b)), lp);
            if ~isempty(r)
                [v, f] = check(r, H, Pn, lidar.covariance(b), c.gate_lidar);
                worst = max(worst, v);
                nbad = nbad + f;
            end
        end
    end
    nis(4) = worst;
    fail(4) = nbad >= c.lidar_beams_fail;
end
if fresh(L.S_DGPS) && dgps.validity
    [r, H] = diffGps_update_model(xn, dgps, gf);
    [nis(5), fail(5)] = check(r, H, Pn, dgps.covariance, c.gate_dgps);
end
m.stamp = stamp(:);

nf = sum(fail);
raw(33) = nf == 1;
raw(34) = nf >= 2;
ok = nis >= 0;
if any(ok)
    metric(35) = (mean(nis(ok)) - 2) * c.frame;
end

d = diag(Pn);
q = xn(10:13);
raw(37) = any(~isfinite(xn)) || any(~isfinite(Pn(:))) || any(d <= 0) || ...
    any(d(1:3) > c.pos_var_max) || any(d(4:6) > c.vel_var_max) || any(d(10:12) > c.att_var_max) || ...
    abs(norm(q) - 1) > c.quat_norm_tol || norm(xn(4:6)) > c.speed_max || norm(xn(7:9)) > c.rate_max;
end

function [v, fail] = check(r, H, P, R, gate)
% normalized NIS (NIS / dof) and whether it fails the 99 % gate
S = H * P * H' + R;
S = 0.5 * (S + S');
n2 = r(:)' * (S \ r(:));
v = n2 / numel(r);
fail = n2 > gate;
end

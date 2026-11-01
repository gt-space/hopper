function [raw, metric, nis] = fdir_nav_monitors(diag, xn, Pn, shadow, c)
% FDIR_NAV_MONITORS  Estimator health checks, catalog rows 33-37.
%
%   diag     the filter's per-sensor update diagnostics this step (16 x 1,
%            gmekf_correct_hopper): (1:8) 1 accepted / 0 gated out / -1 no
%            sample, (9:16) NIS, for baro, gyro, GPS, lidar 1-4, diff GPS
%   xn, Pn   nav estimate (20) and covariance (19 x 19)
%   shadow   the shadow vertical estimator [z; vz; var z; var vz]
%            (shadow_vertical_filter), nav convention
%   c        gates and bounds (gnc_monitor_params)
%
%   raw     37 x 1 logical, rows 33, 34, 36, 37 set (35 is a CUSUM)
%   metric  37 x 1, row 35's CUSUM increment
%   nis     5 x 1 normalized NIS this step (baro gyro gps lidar dgps; -1
%           when that sensor had no update), for logs and tuning
%
% The filter reports each update's NIS (prior, r' S^-1 r) and its own 99 %
% gate decision. FDIR judges the NIS against stricter 99.99 % gates: the
% filter drops single outliers, FDIR looks for a sensor that keeps failing.
%   33  exactly one sensor group fails (lidar: two beams at once)
%   34  two or more fail together (divergence)
%   35  CUSUM of the mean normalized NIS above 2 (poor tuning), per second
%   36  main vs shadow: height or vertical velocity apart by more than
%       c.shadow_sigmas combined sigmas (and the floors c.shadow_dz /
%       c.shadow_dvz)
%   37  covariance or state implausible: not finite, a variance <= 0 or
%       above its bound, quaternion off unit length, speed or rate beyond
%       physical limits
%
%#codegen

raw = false(37, 1); metric = zeros(37, 1); nis = -ones(5, 1);
diag = diag(:); xn = xn(:);
acc = diag(1:8); ni = diag(9:16);
dof = [1 3 6 1 1 1 1 3];
gates = [c.gate_baro c.gate_gyro c.gate_gps c.gate_lidar c.gate_lidar c.gate_lidar c.gate_lidar c.gate_dgps];
upd = acc >= 0 & ni >= 0;
bad = upd(:)' & ni(:)' > gates;

fail = false(5, 1);
groups = {1, 2, 3, 4:7, 8};
for g = 1:5
    k = groups{g};
    if any(upd(k))
        nis(g) = max(ni(k(upd(k))) ./ dof(k(upd(k)))');
    end
end
fail([1 2 3 5]) = bad([1 2 3 8]);
fail(4) = sum(bad(4:7)) >= c.lidar_beams_fail;

nf = sum(fail);
raw(33) = nf == 1;
raw(34) = nf >= 2;
ok = nis >= 0;
if any(ok)
    metric(35) = (mean(nis(ok)) - 2) * c.frame;
end

% main vs shadow vertical estimate
Pd = [Pn(3, 3); Pn(6, 6)];
dz = abs(xn(3) - shadow(1));
dvz = abs(xn(6) - shadow(2));
raw(36) = dz > max(c.shadow_dz, c.shadow_sigmas * sqrt(Pd(1) + shadow(3))) || ...
          dvz > max(c.shadow_dvz, c.shadow_sigmas * sqrt(Pd(2) + shadow(4)));

d = [Pn(1, 1); Pn(2, 2); Pn(3, 3); Pn(4, 4); Pn(5, 5); Pn(6, 6); Pn(10, 10); Pn(11, 11); Pn(12, 12)];
q = xn(10:13);
raw(37) = any(~isfinite(xn)) || any(~isfinite(Pn(:))) || any(diagOf(Pn) <= 0) || ...
    any(d(1:3) > c.pos_var_max) || any(d(4:6) > c.vel_var_max) || any(d(7:9) > c.att_var_max) || ...
    abs(norm(q) - 1) > c.quat_norm_tol || norm(xn(4:6)) > c.speed_max || norm(xn(7:9)) > c.rate_max;
end

function d = diagOf(P)
n = size(P, 1); d = zeros(n, 1);
for i = 1:n, d(i) = P(i, i); end
end

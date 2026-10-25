function result = checkConstraints(result, limits)
%CHECKCONSTRAINTS  Pass/fail for one Monte Carlo result (from mc_main).
%   result = checkConstraints(result)          mission requirements (defaults
%                                              below)
%   result = checkConstraints(result, limits)  override any field of the
%                                              defaults, e.g. the looser CI
%                                              regression limits (mcLimits)
%   landing_vel is the vertical impact speed when the sim stops at ground
%   contact (mc_main), including the drop after engine cutoff.
%   Adds result.flight.landing_offset (m from the pad center at the end of
%   the run) and result.flight.max_tilt_deg (body x-axis from vertical,
%   90 - pitch for this ZYX Euler set), then sets result.status.pass and
%   result.status.failed (names of the violated limits).

L = defaultLimits();
if nargin > 1
    for f = fieldnames(limits)'
        L.(f{1}) = limits.(f{1});
    end
end

ts = result.timeseries;
result.flight.landing_offset = hypot(ts.x(end), ts.y(end));
result.flight.max_tilt_deg   = max(abs(90 - ts.pitch));

v = result.vehicle;
f = result.flight;
checks = {
    'twr_initial',      v.twr_initial    >= L.min_twr
    'dry_mass',         v.dry_mass       >  0
    'wet_mass',         v.wet_mass       >  0
    'max_altitude',     f.max_altitude   >= L.min_altitude
    'max_ascent_vel',   f.max_ascent_vel  <= L.max_ascent_vel
    'max_descent_vel',  f.max_descent_vel <= L.max_descent_vel
    'landing_vel',      f.landing_vel    <= L.max_landing_vel
    'landing_offset',   f.landing_offset <= L.max_landing_offset
    'max_tilt_deg',     f.max_tilt_deg   <= L.max_tilt_deg
    };
ok = [checks{:, 2}];
result.status.pass   = all(ok);
result.status.failed = checks(~ok, 1)';
result.status.limits = L;
end

function L = defaultLimits()
% Mission requirements. (OUT.Vehicle.FinalTWR is minimum throttle over final
% weight, which must be below 1 to descend, so it is not checked against
% min_twr.)
L.min_twr            = 1.0;   % initial, full throttle
L.min_altitude       = 51;    % m
L.max_ascent_vel     = 7;     % m/s
L.max_descent_vel    = 7;     % m/s
L.max_landing_vel    = 0.5;   % m/s, vertical at the end of the run
L.max_landing_offset = 0.5;   % m, pad radius (mc_plots)
L.max_tilt_deg       = 10;    % deg
end

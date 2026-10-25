function L = mcLimits()
%MCLIMITS  Regression limits the CI Monte Carlo enforces (checkConstraints).
%   Today's vehicle plus margin, so CI catches a change that makes flights
%   worse. These are looser than the mission requirements (checkConstraints
%   defaults), which runMonteCarlo reports separately as meets_requirements.
%   Measured 2026-10-24, seeded 10-run MC: ascent 8.13, descent 7.35,
%   impact 1.03-1.04 m/s, offset 0.73-2.87 m, tilt 0.83-2.09 deg.
%   Tighten these as the vehicle improves; note any change in the PR.

L.min_altitude       = 51;    % m
L.max_ascent_vel     = 8.5;   % m/s (requirement 7; Trajectory2 peaks at 8.125)
L.max_descent_vel    = 7.7;   % m/s (requirement 7)
L.max_landing_vel    = 1.2;   % m/s impact (requirement 0.5)
L.max_landing_offset = 3.5;   % m (requirement 0.5, pad radius)
L.max_tilt_deg       = 5;     % deg
end

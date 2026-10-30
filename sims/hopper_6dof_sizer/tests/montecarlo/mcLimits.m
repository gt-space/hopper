function L = mcLimits()
%MCLIMITS  Regression limits the CI Monte Carlo enforces (checkConstraints).
%   Today's vehicle plus margin, so CI catches a change that makes flights
%   worse. These are looser than the mission requirements (MISSION CONSTRAINTS
%   in inputs/mission_inputs.m), which runMonteCarlo reports separately as
%   meets_requirements.
%   Measured 2026-10-24, seeded 10-run MC: ascent 8.13, descent 7.35,
%   impact 1.03-1.04 m/s, offset 0.73-2.87 m, tilt 0.83-2.09 deg.
%   From the 2026-10-29 seeded 10-run MC (outputs/montecarlo/mc_ci.mat):
%   touchdown lateral speed 0.003-0.032 m/s, touchdown tilt 0.24-1.15 deg.
%   Tighten these as the vehicle improves; note any change in the PR.

L.min_altitude       = 51;    % m
L.max_ascent_vel     = 8.5;   % m/s (requirement 10; Trajectory2 peaks at 8.125)
L.max_descent_vel    = 7.7;   % m/s (requirement 10)
L.max_landing_vel    = 1.2;   % m/s impact (requirement 1.5)
L.max_lateral_vel    = 0.1;   % m/s at touchdown (requirement 0.5)
L.max_landing_offset = 3.5;   % m (requirement 1)
L.max_tilt_deg       = 5;     % deg (requirement 10)
L.max_touchdown_tilt_deg = 2; % deg (requirement 2)
end

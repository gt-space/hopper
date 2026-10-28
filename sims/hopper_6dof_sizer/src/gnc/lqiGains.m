function [Klqi, p] = lqiGains()
%LQIGAINS  Gains for the CONTROLLER=lqi choice (LQI-DEV's hover regulator).
%   [Klqi, p] = lqiGains() returns the 4x15 LQI gain (12 states + 3 position
%   integrators) and the trim attitude as MRPs, the two variables the
%   harness's LQI choice reads.
%
%   Same computation as LQIDevR1 (LQI-DEV branch, "half finished"), as a
%   function: no clear all, and the gain comes back as Klqi, the name the
%   v4 model reads (LQIDevR1 names it KlqI). LQI-DEV linearizes a unit
%   plant (m = 1 kg, I = eye(3), hover thrust 10 N), so these gains are
%   not tuned for the hopper or the CTV. Retune here, not in the model.

Qdiags  = [0 0 0 10 10 10 1e2 1e2 1e2 1e6 1e6 1e6];
QIdiags = [10 10 500];
QI = diag([Qdiags, QIdiags]);
R  = diag([100, 1e8, 1e8, 1e4]);     % T, delta_y, delta_z, rcs

params = struct('g', 9.81, 'd', 0.2, 'r_m', 0.125, 'n_rcs', 1, 'I', eye(3));
m = 1;

q = [sqrt(2)/2; 0; sqrt(2)/2; 0];    % trim: body x up
p = qToMRP(q);
x = [zeros(9, 1); p];
u = [10; 0; 0; 0];

[A, B] = HopperLinearization_lqi(x, u, params, m);
[Atilde, Btilde] = extendedDynamics(A, B);
Klqi = lqr(Atilde, Btilde, QI, R, zeros(15, 4));
end

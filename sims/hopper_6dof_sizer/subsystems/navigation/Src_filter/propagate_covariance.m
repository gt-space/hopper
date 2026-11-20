%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% propagate_covariance.m
% 
% Function to propagate the covariance over a time step.
% 
% Author: David Meisinger
% Date: 10/11/26
% 
% Inputs:
%   Pn: n covariance (19x19)
%   F05: discrete error-state transition matrix (19x19)
%   G05: discrete process-noise Jacobian (19x13)
%   Qc: per-sample IMU covariance and continuous-time bias random-walk PSDs
%   dt: timestep (s), used for input validation
% 
% Outputs:
%   Pnp1: n+1 covariance (19x19)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pnp1 = propagate_covariance(Pn, F05, G05, Qc, dt)

    % Size checks
    assert(dt >= 0);
    assert(all(size(Pn) == [19,19]));
    assert(all(size(F05) == [19,19]));
    assert(all(size(G05) == [19,13]));
    assert(all(size(Qc) == [13,13]));

    % The IMU noise is sampled once per propagation; bias random walks are
    % scaled by sqrt(dt) in G05.
    Pnp1 = F05 * Pn * F05' + G05 * Qc * G05';
    Pnp1 = 0.5 * (Pnp1 + Pnp1');

end
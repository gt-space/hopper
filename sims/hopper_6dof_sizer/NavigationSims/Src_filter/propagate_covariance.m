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
%   F05: state Jacobian for n+0.5 time
%   G05: noise Jacobian for n+0.5 time
%   Qc: process noise intensity matrix
%   dt: timestep (s)
% 
% Outputs:
%   Pnp1: n+1 covariance (19x19)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pnp1 = propagate_covariance(Pn, F05, G05, Qc, dt)

    stateDim = size(Pn,1);
    
    % Size checks
    assert(dt >= 0);
    assert(all(size(Pn) == [19,19]));
    assert(all(size(F05) == [19,19]));
    assert(all(size(G05) == [19,13]));
    assert(all(size(Qc) == [13,13]));

    % Second-Order-Taylor-Expansion of STM
    STM = eye(stateDim) + F05 * dt + ...
        0.5 * (F05 * F05) * dt^2;

    %{
        Half-Step STM for evaluating process noise
        integration at the midpoint.
    %}
    halfStep = dt / 2;
    STMhalf = eye(stateDim) + F05 * halfStep + 0.5 * (F05 * F05) * halfStep^2;

    %{
        Projects continuos-time process noise intensity
        to state-space through noise Jacobian.
    %}
    processNoiseIntensity = G05 * Qc * G05';
    
    % Midpoint approximation of discrete process-noise covariance
    discreteProcessNoiseCovariance = STMhalf * processNoiseIntensity * STMhalf' * dt;

    % Propagation and Symmetrization
    Pnp1 = STM * Pn * STM' + discreteProcessNoiseCovariance;
    Pnp1 = 0.5 * (Pnp1 + Pnp1');

end
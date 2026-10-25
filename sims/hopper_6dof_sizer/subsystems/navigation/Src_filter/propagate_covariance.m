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
    
    % Estimate STM to second order
    STM = eye(19) + F05*dt + 1/2*F05^2*dt^2;

    % Estimate process noise covariance (first order overestimate)
    P_ww = STM*G05*Qc*G05'*STM'*dt;
    
    % Propagate covariance
    Pnp1 = STM*Pn*STM' + P_ww;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update.m
% 
% Generic function to update a state estimate using a given sensor model
% and measurement.  
% 
% Author: David Meisinger
% Date: 10/18/26
% 
% Inputs:
%   x: state vector
%   P: covariance
%   model: function handle for sensor model
%   meas: measurement struct for the sensor
%   params: Any params needed for the sensor model
% 
% Outputs:
%   x_c: corrected state vector
%   P_c: corrected covariance
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_c, P_c] = update(x, P, model, meas, params)

    % Get the measurement innovation (estimate - meas) and Jacobian column
    [r, H] = model(x, meas, params);
    R = meas.R;
    S = H*P*H' + R;

    % chi-square innovation gate
    if r'*(S\r) > gate(numel(r))
        x_c = x; P_c = P;
        return;
    end
    
    K  = P*H'/S;
    dx = K*r;                              % 19x1 error estimate (est - truth)
    I_KH = eye(19) - K*H;
    P_c = I_KH*P*I_KH' + K*R*K';            % Joseph form
    P_c = 0.5*(P_c+P_c');
    x_c = inject(x, dx);                     % corrections subtract, since est - truth
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmekf_predict.m
% 
% Function to generate the propagated state and covariance estimates used
% in the GMEKF filter.  
% 
% Author: David Meisinger
% Date: 10/18/26
% 
% Inputs:
%   xn: n state vector (20x1)
%   Pn: n covariance (19x19)
%   accel_meas: IMU acceleration measured (3x1) (m/s^2) (body frame)
%   angvel_meas: IMU angular velocity measured (3x1) (deg/s) (body frame)
%   validity: if IMU measurement is valid
%   J: Inertia tensor, xyz aligning with body frame (kg*m^2)
%   g: acceleration of gravity (scalar)
%   dt: time delta to propagate (s)
%   Qc: process noise intensity matrix (13x13)
% 
% Outputs:
%   xnp1: n+1 state (20x1)
%   Pnp1: n+1 covariance (19x19)
%   aux: struct with useful shared variables for the filter
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnp1, Pnp1, aux] = gmekf_predict(xn, Pn, accel_meas, angvel_meas, validity, J, g, dt, Qc)
    
    % Generate some shared quantities
    % NOTE: aux is built in one go at each stage, since code generation
    % cannot add fields to a struct after it has been used
    g_vec = [0; 0; g];
    A_IB = quaternion_to_rotation(xn(10:13));
    auxProp = struct('A_IB', A_IB);

    % Propagate the state
    xnp1 = propagate_state(xn, auxProp, accel_meas, angvel_meas, validity, J, g, dt);
    xnp05 = propagate_state(xn, auxProp, accel_meas, angvel_meas, validity, J, g, dt*0.5);
    % ^xnp05 is the propagation to time step 0.5 which allows better STM
    % estimation

    % Calculate rotation matrix at the one-half time step
    A_IB05 = quaternion_to_rotation(xnp05(10:13));

    % fill aux params
    if validity
        a_I = A_IB05*(accel_meas - xnp05(14:16));
    else
        a_I = g_vec; % Assume vehicle is unaccelerated
    end
    w_I = A_IB05*xnp05(7:9);

    % Get the error state Jacobians
    auxJac = struct('A_IB', A_IB, 'A_IB05', A_IB05, 'a_I', a_I, 'w_I', w_I);
    [F05, G05] = error_state_jacobian(xnp05, auxJac, J);
    aux = struct('A_IB', A_IB, 'A_IB05', A_IB05, 'a_I', a_I, 'w_I', w_I, 'F', F05, 'G', G05);

    % Propagate the covariance
    Pnp1 = propagate_covariance(Pn, F05, G05, Qc, dt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lidar_system_model.m
% 
% Function to simulate the measurement of four LIDAR sensors given the
% vehicle state and sensor frames.
% Depends on filter_sim_initializer setting the sensor frames.
% 
% Author: David Meisinger
% Date: 10/17/26
% 
% Inputs:
%   x: true state vector (20x1)
% 
% Outputs:
%   d_vec: LIDAR sensor distance estimates (4x1) (m)
%   timestamp: time elapsed (s)
%   validity: boolean vector (true if valid measurement) (4x1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d_vec, timestamp, validity] = lidar_system_model(x, time, lidar_frames)

    % Load the LIDAR frame variables
    A1 = lidar_frames.A1;
    A2 = lidar_frames.A2;
    A3 = lidar_frames.A3;
    A4 = lidar_frames.A4;
    r1 = lidar_frames.r1;
    r2 = lidar_frames.r2;
    r3 = lidar_frames.r3;
    r4 = lidar_frames.r4;

    % Get the height and attitude of the vehicle
    h = x(3);
    q_BI = x(10:13);
    q_IB = [q_BI(1); -q_BI(2:4)];

    % Compute the rotation matrix 3rd row (n) from q components
    q0 = q_IB(1); q1 = q_IB(2); q2 = q_IB(3); q3 = q_IB(4);
    n = [2*(q3*q1+q2*q0) 2*(q3*q2-q1*q0) -q1^1-q2^2+q3^2+q0^2];

    % Compute all the LIDAR true distances
    d_vec = [-(n*r1 + h)/(n*A1(:,1));
             -(n*r2 + h)/(n*A2(:,1));
             -(n*r3 + h)/(n*A3(:,1));
             -(n*r4 + h)/(n*A4(:,1))];

    % Apply sensor model with noise
    [d_vec, validity] = lidar_model_base(d_vec);

    % Set timestamp (assume no delay for now)
    timestamp = time;

end
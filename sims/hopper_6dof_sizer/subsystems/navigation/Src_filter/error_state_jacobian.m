%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error_state_jacobian.m
%
% Discrete error-state transition and noise Jacobian for one IMU sample.
%
% Inputs:
%   x: estimated state at the propagation midpoint (20x1)
%   aux: struct containing A_IB, A_IB05, a_I, and w_I
%   dt: timestep (s)
%   validity: whether the IMU sample is valid
%   J: body-frame inertia tensor
%
% Outputs:
%   F: discrete error-state transition matrix (19x19)
%   G: discrete process-noise Jacobian (19x13)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, G] = error_state_jacobian(x, aux, dt, validity, J)

    F = eye(19);
    G = zeros(19,13);

    C = aux.A_IB05;
    a_I = aux.a_I;
    w_I = aux.w_I;

    F(1:3, 4:6) = dt*eye(3);

    if validity
        F(1:3, 10:12) = -0.5*skew(a_I)*dt^2;
        F(1:3, 13:15) = -0.5*eye(3)*dt^2;
        F(4:6, 10:12) = -skew(a_I)*dt;
        F(4:6, 13:15) = -eye(3)*dt;

        % The measured rate replaces the previous rate and depends directly
        % on the body-frame gyro bias; its error is inertial-aligned.
        F(7:9, 7:9) = zeros(3);
        F(7:9, 16:18) = -aux.A_IB';

        % Attitude uses the mean of the old rate and the bias-corrected
        % sample rate.
        F(10:12, 7:9) = 0.5*C*dt;
        F(10:12, 16:18) = -0.5*eye(3)*dt;

        % Accelerometer and gyro sample noise; Qc contains per-sample
        % measurement covariance in its first six channels.
        G(1:3, 1:3) = 0.5*C*dt^2;
        G(4:6, 1:3) = C*dt;
        G(7:9, 4:6) = eye(3);
        G(10:12, 4:6) = 0.5*C*dt;
    else
        angvel = x(7:9);
        M = J \ (skew(J*angvel) - skew(angvel)*J);
        F(7:9, 7:9) = eye(3) + M*dt;
        F(10:12, 7:9) = C*(eye(3) + 0.5*M*dt)*dt;
    end

    % Bias errors are represented in inertial coordinates while the nominal
    % biases are body-frame quantities.
    F(13:15, 13:15) = eye(3) + skew(w_I)*dt;
    F(16:18, 16:18) = eye(3) + skew(w_I)*dt;

    % Bias random walks and barometer-bias random walk are continuous-time
    % spectral densities, scaled here to one discrete sample.
    G(13:15, 7:9) = sqrt(dt)*aux.A_IB05;
    G(16:18, 10:12) = sqrt(dt)*aux.A_IB05;
    G(19, 13) = sqrt(dt);

    assert(isequal(size(x), [20 1]));
end

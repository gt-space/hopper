% This script must be executed prior to running the Simulink model
% "Hopper_6dof_nonlinear_model_openLoop" in order to initialize the
% integrator states and set the model’s initial conditions.

% State vector initial conditions:
x0 = [ 0; 0; 0;         % Position: X, Y, Z (NED). Example: Z = -50 = 50 m altitude.
          0; 0; 0;          % Linear velocities: u, v, w (body frame)
          0; 0; 0;          % Angular rates: p, q, r (body frame)
        0.7071; 0; 0.7071; 0];     % Quaternion attitude: [q0 q1 q2 q3].
                                % q = [1 0 0 0] aligns the body frame with the NED frame.
        
                                 % To align the body X-axis with the NED Z-axis (vehicle nose
                                 % pointing downward in NED), use the quaternion:
                                 %     q = [0.7071; 0; 0.7071; 0];

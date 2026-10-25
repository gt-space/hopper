function resetJacobian = compute_error_reset_jacobian(errorCorrection)
% (Yuma's implementation | 10/24/26)
% Computes the local error-reset Jacobian for covariance update
%
% Inputs:
%
%   errorCorrection - Estimated error correction state vector (19x1)
%
% Outputs:
%
%   resetJacobian   - Local covariance reset matrix (19x19)  

    assert(numel(errorCorrection) == 19);

    attitudeError = errorCorrection(10:12);
    accelBiasError = errorCorrection(13:15);
    gyroBiasError = errorCorrection(16:18);

    resetJacobian = eye(19);
    
    % First-order attitude-error reset
    resetJacobian(10:12, 10:12) = eye(3) - 0.5 * skew(attitudeError);

    % Coupling attitude and inertial-aligned bias errors
    resetJacobian(13:15,10:12) = -skew(accelBiasError);
    resetJacobian(16:18,10:12) = -skew(gyroBiasError);

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute_error_reset_jacobian.m
% Computes the local error reset Jacobian matrix J for the covariance update.
% 
% Inputs:
%   delta_x_hat: Estimated error state vector correction
% 
% Outputs:
%   J: Local covariance reset matrix (19x19 or matching error state dimension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
function J = compute_error_reset_jacobian(delta_x_hat)
delta_theta = delta_x_hat(1:3);
% Skew-symmetric matrix of attitude error
theta_skew = [0, -delta_theta(3), delta_theta(2); ...
    delta_theta(3), 0, -delta_theta(1); ...
    -delta_theta(2), delta_theta(1), 0];

% First-order reset Jacobian for MEKF: J_M = diag(I_3 + 0.5*[theta]_x, I_3, I_3, I_3, I_3, 1)[cite: 4]
dim = length(delta_x_hat);
J = eye(dim);
% Apply attitude error coupling to the first 3x3 block
J(1:3, 1:3) = eye(3) + 0.5 * theta_skew;
% Remaining blocks are identity mappings for velocity, position, and biases to first order[cite: 4]
end

%}

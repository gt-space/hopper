%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state_injection.m
% Injects the estimated error correction vector into the nominal state.
% 
% Inputs:
%   m_minus: Prior nominal state vector
%   delta_x_hat: Estimated error state vector (correction)
% 
% Outputs:
%   m_plus: Updated nominal state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function m_plus = state_injection(m_minus, delta_x_hat)
% Extract error states from delta_x_hat
delta_theta = delta_x_hat(1:3); % Attitude error
delta_v     = delta_x_hat(4:6); % Velocity error
delta_chi   = delta_x_hat(7:9); % Position error
delta_beta  = delta_x_hat(10:12);% Gyro bias error
delta_gamma = delta_x_hat(13:15);% Accel bias error
delta_lambda= delta_x_hat(16);   % Baro bias error (or adjust indices to match your state layout)

% Extract nominal states (Adjust indices to match your specific state vector layout)
% Assuming nominal state order: [Position(3), Velocity(3), Quaternion(4), GyroBias(3), AccelBias(3), BaroBias(1)]
m_X   = m_minus(1:3);
m_V   = m_minus(4:6);
m_q   = m_minus(7:10);
m_bg  = m_minus(11:13);
m_ba  = m_minus(14:16);
m_bh  = m_minus(17);

% --- Apply Corrections (MEKF Injection Rules) ---
% 1. Position update: m_X^+ = m_X^- - delta_chi[cite: 4]
m_X_plus = m_X - delta_chi;

% 2. Velocity update: m_V^+ = m_V^- - delta_v[cite: 4]
m_V_plus = m_V - delta_v;

% 3. Attitude update via quaternion exponential map: m_q^+ = m_q^- \otimes Exp_q(-\delta\theta)[cite: 4]
theta_norm = norm(delta_theta);
if theta_norm < 1e-8
    dq = [0; 0; 0; 1];
else
    half_theta = -0.5 * delta_theta;
    dq = [sin(norm(half_theta)) * (half_theta / norm(half_theta)); cos(norm(half_theta))];
end
m_q_plus = quatmultiply(m_q', dq')'; % Standard quaternion multiplication
m_q_plus = m_q_plus / norm(m_q_plus); % Normalize

% 4. Bias updates[cite: 4]
m_bg_plus  = m_bg - delta_beta;
m_ba_plus  = m_ba - delta_gamma;
m_bh_plus  = m_bh - delta_lambda;

% Recombine into m_plus state vector
m_plus = [m_X_plus; m_V_plus; m_q_plus; m_bg_plus; m_ba_plus; m_bh_plus];
end
%}


function nominalStatePosterior = state_injection(nominalStatePrior, errorCorrection)
% (Yuma's implementation | 10/23/26)
% Measurement update for state estimate
%
% Inputs:
% 
%   nominalStatePrior       - Prior state estimate (20x1)
%   errorCorrection         - Estimated error correction state vector (19x1)
% 
% Outputs:
%
%   nominalStatePosterior   - Posterior state estimate (20x1)
    
    % Extract nominal state
    positionPrior        = nominalStatePrior(1:3);
    velocityPrior        = nominalStatePrior(4:6);
    angularVelocityPrior = nominalStatePrior(7:9);
    quaternionPrior      = nominalStatePrior(10:13);
    accelBiasPrior       = nominalStatePrior(14:16);
    gyroBiasPrior        = nominalStatePrior(17:19);
    baroBiasPrior        = nominalStatePrior(20);
    
    % Extract estimated error components
    positionError        = errorCorrection(1:3);
    velocityError        = errorCorrection(4:6);
    angularVelocityError = errorCorrection(7:9);
    attitudeError        = errorCorrection(10:12);
    accelBiasError       = errorCorrection(13:15);
    gyroBiasError        = errorCorrection(16:18);
    baroBiasError        = errorCorrection(19);

    positionPosterior = positionPrior - positionError;
    velocityPosterior = velocityPrior - velocityError;
    baroBiasPosterior = baroBiasPrior - baroBiasError;
    
    % Remaining injection equations to be derived
end
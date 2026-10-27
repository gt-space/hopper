function nominalStatePosterior = state_injection(nominalStatePrior, errorCorrection)
% (Yuma's implementation | 10/24/26)
% Injects the estimated error correction vector into the nominal state
%
% Inputs:
%
%   nominalStatePrior       - Prior state estimate (20x1)
%   errorCorrection         - Estimated error correction state vector (19x1)
% 
% Outputs:
%
%   nominalStatePosterior   - Posterior state estimate (20x1)
%
% NOTES:
%   Position and velocity are expressed in the inertial NED frame
%   Nominal angular velocity and IMU biases are body-frame quantities
%   GMEKF accelerometer and gyro bias errors are inertial-aligned
    
    % Extract nominal state components
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

    % Additive corrections for inertial position and velocity errors
    positionPosterior = positionPrior - positionError;
    velocityPosterior = velocityPrior - velocityError;
    
    %{
        Correct the quaternion using the attitude error small-angle approximation.
        For the full mathematical derivation, see:
        https://app.notion.com/p/GMEKF-Navigation-Filter-Conventions-Formulas-Used-40327990d78e803d95c4ec19b3e22e12

        FOR FUTURE: May consider finite rotation formula with exponential-maps attitude correction
    %}
    quaternionPosterior = quaternionPrior - 0.5 * Psi_helper(quaternionPrior) * attitudeError;

    % Re-normalize the quaternion
    quaternionPosterior = quaternionPosterior / norm(quaternionPosterior);


    rotationBodyToInertialPosterior = quaternion_to_rotation(quaternionPosterior);

    % Frame corrections
    accelBiasPosterior = accelBiasPrior - rotationBodyToInertialPosterior' * accelBiasError;
    gyroBiasPosterior = gyroBiasPrior - rotationBodyToInertialPosterior' * gyroBiasError;

    baroBiasPosterior = baroBiasPrior - baroBiasError; % correct additive barometer bias
    angularVelocityPosterior = angularVelocityPrior - angularVelocityError; % correct additive angular velocity error


    % Assemble the state
    nominalStatePosterior = ...
        [   positionPosterior;
            velocityPosterior;
            angularVelocityPosterior;
            quaternionPosterior;
            accelBiasPosterior;
            gyroBiasPosterior;
            baroBiasPosterior   ];

end
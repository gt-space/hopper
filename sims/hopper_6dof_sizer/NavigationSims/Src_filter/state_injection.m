%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state_injection.m
% Injects the estimated error correction vector into the nominal state.
% 
% Author: Yuma Nagayoshi
% Date: 10/23/26
% 
% Inputs:
%   nominalStatePrior       - Prior state estimate (20x1)
%   errorCorrection         - Estimated error correction state vector (19x1)
% 
% Outputs:
%   nominalStatePosterior   - Posterior state estimate (20x1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nominalStatePosterior = state_injection(nominalStatePrior, errorCorrection)
    
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

    % Make the additive corrections
    positionPosterior = positionPrior - positionError;
    velocityPosterior = velocityPrior - velocityError;
    angularVelocityPosterior = angularVelocityPrior + angularVelocityError;
    accelBiasPosterior = accelBiasPrior + accelBiasError;
    gyroBiasPosterior = gyroBiasPrior + gyroBiasError;
    baroBiasPosterior = baroBiasPrior - baroBiasError;
    
    % Correct the quaternion using the attitude error, which is a small
    % angle rotation vector
    w = quaternionPrior(1); x = quaternionPrior(2);
    y = quaternionPrior(3); z = quaternionPrior(4);
    Chi = [-x -y -z;
            w -z  y;
            z  w -x;
           -y  x  w];
    quaternionPosterior = quaternionPrior + 0.5*Chi*attitudeError;

    % Renormalize the quaternion
    quaternionPosterior = quaternionPosterior/norm(qauternionPosterior);

    % Assemble the state
    nominalStatePosterior = [positionPosterior;
                             velocityPosterior;
                             angularVelocityPosterior;
                             quaternionPosterior;
                             accelBiasPosterior;
                             gyroBiasPosterior;
                             baroBiasPosterior];
end
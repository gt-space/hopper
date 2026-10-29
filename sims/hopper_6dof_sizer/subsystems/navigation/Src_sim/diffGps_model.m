% diffGps_model.m
%
% Function to calculate a noisy differential GPS measurement operating on
% the true vehicle state. The measurement is the baseline vector from
% antenna 1 to antenna 2, expressed in the inertial NED frame.
%
% Coordinate conventions:
%   Baseline measurement: inertial NED
%   Quaternion: inertial-to-body, Hamilton, scalar-first
%   Antenna positions: body frame, relative to COM
%
% Noise model:
%   Carrier-phase relative positioning error is roughly constant in
%   metres, so the angular error scales with 1/(baseline length). Noise is
%   applied per NED axis; heading accuracy follows from the baseline.
%
% Inputs:
%   trueState       - True vehicle state (20x1)
%   simtime         - Simulation time (s)
%   gpsFrameBus     - Struct with antenna placement from COM
%
% Outputs:
%   diffGpsBus      - struct defined in the Navigation Dictionary
function diffGpsBus = diffGps_model(trueState, simtime, gpsFrameBus)

    % Load in the diffGps_params struct
    persistent diffGps_params;
    if isempty(diffGps_params)
        diffGps_params.noise_std = [0.007; 0.007; 0.012];          % Baseline jitter per NED axis (m)
        diffGps_params.covariance = diag([0.007; 0.007; 0.012].^2); % Filter measurement covariance
        diffGps_params.packet_invalid_prob = 0.0;                  % Set to 0 for now
    end

    timestamp = simtime;

    quaternionTrue = trueState(10:13);

    % Body-to-inertial NED rotation
    C_IB = quaternion_to_DCM(quaternionTrue);

    % Antenna 1 -> antenna 2 baseline in body coordinates
    baselineBody = gpsFrameBus.r2(:) - gpsFrameBus.r1(:);

    % Ideal baseline in NED
    baselineInertial = C_IB * baselineBody;

    % Measurement noise
    baselineNoise = diffGps_params.noise_std(:) .* randn(3,1);

    % Simulated measurement
    basisVec = baselineInertial + baselineNoise;
    validity = rand() >= diffGps_params.packet_invalid_prob;

    % Load the output struct
    diffGpsBus.basisVec = basisVec;
    diffGpsBus.covariance = diffGps_params.covariance;
    diffGpsBus.timestamp = timestamp;
    diffGpsBus.validity = validity;

end

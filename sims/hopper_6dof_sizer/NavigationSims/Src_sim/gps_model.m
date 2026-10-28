% GPS antenna position and velocity simulation.
%
% Coordinate conventions:
%   Position, velocity: inertial NED
%   Angular velocity: body frame
%   Quaternion: inertial-to-body, Hamilton, scalar-first
%   Antenna lever arm: body frame
%
% Inputs:
%   trueState       - True vehicle state (20x1)
%   simtime         - Simulation time (s)
%   gpsFrameBus     - Struct with antenna placement from COM
%
% Outputs:
%   gpsBus          - struct defined in the Navigation Dictionary
function gpsBus = gps_model(trueState, simtime, gpsFrameBus)

    % Load in the gps_params struct
    persistent gps_params;
    if isempty(gps_params)
        gps_params.r1 = [0.25; 0.0; -0.2];             % Antenna 1 vector from COM (3x1) (m)
        gps_params.r2 = [-0.25; 0.0; -0.2];            % Antenna 2 vector from COM (3x1) (m)
        gps_params.bias = [0; 0; 0];                   % Constant position bias (m)
        gps_params.multipath = [0.002; 0.001; -0.003]; % Pad bounce error (b_mp,k) (m)[cite: 1]
        gps_params.pos_noise_std = [0.02; 0.02; 0.05]; % Position jitter standard deviation (eta_p) (m)[cite: 1]
        gps_params.vel_noise_std = [0.05, 0.05, 0.01]; % Velocity jitter standard deviation
        gps_params.pos_covariance = diag([0.02^2; 0.02^2; 0.05^2]); % Filter measurement covariance
        gps_params.vel_covariance = diag([0.05^2; 0.05^2; 0.1^2]); % Filter measurement covariance
        gps_params.packet_invalid_prob = 0.0;          % Set to 0 for now
        % Provisional GPS velocity noise standard deviations (m/s).
        % Replace with receiver specifications.
        gps_params.vel_noise_std = [0.1; 0.1; 0.1];
    end
    
    timestamp = simtime;
    
    positionTrue = trueState(1:3);
    velocityTrue = trueState(4:6);
    angularVelocityTrue = trueState(7:9);
    quaternionTrue = trueState(10:13);
    
    % Body-to-inertial NED rotation
    C_IB = quaternion_to_DCM(quaternionTrue);
    
    % Antenna lever arm in body coordinates (point the receiver reports)
    r_B = gps_lever_arm(gpsFrameBus);
    
    % Ideal antenna position and velocity
    antennaPosition = positionTrue + C_IB * r_B;
    antennaVelocity = velocityTrue + C_IB * cross(angularVelocityTrue, r_B);
    
    % Measurement noise
    positionNoise = gps_params.pos_noise_std(:) .* randn(3,1);
    velocityNoise = gps_params.vel_noise_std(:) .* randn(3,1);
    
    % Simulated measurements
    position = antennaPosition + gps_params.bias(:) + gps_params.multipath(:) + positionNoise;
    velocity = antennaVelocity + velocityNoise;
    validity = rand() >= gps_params.packet_invalid_prob;

    % Load the output struct
    gpsBus.pos = position;
    gpsBus.vel = velocity;
    gpsBus.pcovariance = gps_params.pos_covariance;
    gpsBus.vcovariance = gps_params.vel_covariance;
    gpsBus.timestamp = timestamp;
    gpsBus.validity = validity;

end
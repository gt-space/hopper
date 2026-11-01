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

    % Simulink gets touchy about dims
    trueState = trueState(1:20,1);

    % Load in the gps_params struct
    persistent gps_params;
    if isempty(gps_params)
        P = nav_sensor_params();   % defaults there; tests can override
        gps_params = P.gps;
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
    correction = cross(angularVelocityTrue, r_B);
    antennaVelocity = velocityTrue + C_IB * correction(1:3,1); % Simulink gets touchy about dims
    
    % Measurement noise
    positionNoise = gps_params.pos_noise_std(:) .* randn(3,1);
    velocityNoise = gps_params.vel_noise_std(:) .* randn(3,1);
    
    % Simulated measurements
    position = antennaPosition + gps_params.bias(:) + gps_params.multipath(:) + positionNoise;
    velocity = antennaVelocity + velocityNoise;
    validity = rand() >= gps_params.packet_invalid_prob;

    % Load the output struct
    gpsBus.pos = position(1:3,1); % Simulink gets touchy about dims
    gpsBus.pcovariance = gps_params.pos_covariance;
    gpsBus.vel = velocity(1:3,1);
    gpsBus.vcovariance = gps_params.vel_covariance;
    gpsBus.timestamp = timestamp;
    gpsBus.validity = validity;

end
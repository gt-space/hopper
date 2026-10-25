function [position, velocity, timestamp, validity] = gps_model(trueState, simtime)
% GPS antenna position and velocity simulation.
%
% Coordinate conventions:
%   Position, velocity: inertial NED
%   Angular velocity: body frame
%   Quaternion: inertial-to-body, Hamilton, scalar-first
%   Antenna lever arm: body frame

    persistent gps_params;

    if isempty(gps_params)
        try
            gps_params = evalin('base', 'gps_params');
        catch
            gps_params.r_b = [0; 0; 0];
            gps_params.bias = [0; 0; 0];
            gps_params.multipath = [0; 0; 0];
            gps_params.pos_noise_std = [0.02; 0.02; 0.05];
            gps_params.vel_noise_std = [0.1; 0.1; 0.1];
            gps_params.packet_invalid_prob = 0;
        end
    end
    
    timestamp = simtime;
    
    positionTrue = trueState(1:3);
    velocityTrue = trueState(4:6);
    angularVelocityTrue = trueState(7:9);
    quaternionTrue = trueState(10:13);
    
    % Body-to-inertial NED rotation
    C_IB = quaternion_to_DCM(quaternionTrue)';
    
    % Antenna lever arm in body coordinates
    r_B = gps_params.r_b(:);
    
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

end
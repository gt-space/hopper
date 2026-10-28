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

    % Build an empty output struct from the dictionary definition
    persistent gpsBusData;
    if isempty(gpsBusData)
        dictionary = Simulink.data.dictionary.open('NavigationDictionary.sldd');
        gpsBusData = Simulink.Bus.createMATLABStruct('gpsBus',[],1,dictionary);
    end

    % Load in the gps_params struct
    % persistent gps_params;
    % if isempty(gps_params)
    %     try
    %         gps_params = evalin('base', 'gps_params');
    %     catch
    %         gps_params.bias = [0; 0; 0];
    %         gps_params.multipath = [0; 0; 0];
    %         gps_params.pos_noise_std = [0.02; 0.02; 0.05];
    %         gps_params.vel_noise_std = [0.1; 0.1; 0.1];
    %         gps_params.pos_covariance = diag([0.02^2; 0.02^2; 0.05^2]); % Filter measurement covariance
    %         gps_params.vel_covariance = diag([0.05^2; 0.05^2; 0.1^2]); % Filter measurement covariance
    %         gps_params.packet_invalid_prob = 0;
    %     end
    % end
    
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
    gpsBusData.pos = position;
    gpsBusData.vel = velocity;
    gpsBusData.pcovariance = gps_params.pos_covariance;
    gpsBusData.vcovariance = gps_params.vel_covariance;
    gpsBusData.timestamp = timestamp;
    gpsBusData.validity = validity;
    gpsBus = gpsBusData;

end
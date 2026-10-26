% IMU measurement simulation.
%
% Inputs:
%   x          - True vehicle state (20x1)
%   a          - True inertial acceleration in NED (3x1)
%   accelBias  - True accelerometer bias, body frame (3x1)
%   gyroBias   - True gyroscope bias, body frame (3x1)
%   simtime    - Simulation time (s)
%
% Outputs:
%   accelBus   - struct defined in the Navigation Dictionary
%   angvelBus  - struct defined in the Navigation Dictionary
function [accelBus, angvelBus] = imu_model(...
    trueState, accelerationTrue, accelerometerBiasTrue, gyroBiasTrue, simtime)
    
    % Build empty output structs
    persistent accelBusData angvelBusData;
    if isempty(accelBusData) || isempty(angvelBusData)
        dictionary = Simulink.data.dictionary.open('NavigationDictionary.sldd');
        accelBusData = Simulink.Bus.createMATLABStruct('accelBus',[],1,dictionary);
        angvelBusData = Simulink.Bus.createMATLABStruct('angvelBus',[],1,dictionary);
    end
    
    % Retrieve simulation parameters
    persistent imu_params;
    if isempty(imu_params)
        try
            imu_params = evalin('base', 'imu_params');
        catch
            imu_params.accel_noise_std = [0; 0; 0];
            imu_params.gyro_noise_std = [0; 0; 0];
            imu_params.gravity_scalar = 9.81;
            imu_params.packet_invalid_prob = 0.0;
            imu_params.acovariance = diag([0.05; 0.05; 0.05].^2);
            imu_params.gcovariance = diag([0.005; 0.005; 0.005].^2);
        end
    end
    
    timestamp = simtime;

    angularVelocityTrue = trueState(7:9);
    quaternionTrue = trueState(10:13);

    rotationInertialToBody = quaternion_to_rotation(quaternionTrue);
    
    % Gravity in inertial NED coordinates
    gravityInertial = [0; 0; imu_params.gravity_scalar];
    
    % True specific force in body coordinates
    specificForceBodyTrue = rotationInertialToBody * (accelerationTrue(:) - gravityInertial);
    

    accelerometerNoise = imu_params.accel_noise_std .* randn(3, 1);
    gyroNoise = imu_params.gyro_noise_std .* randn(3,1);
    
    accelerationSimulated = specificForceBodyTrue + accelerometerBiasTrue(:) + accelerometerNoise;
    angularVelocitySimulated = angularVelocityTrue + gyroBiasTrue(:) + gyroNoise;
    
    
    % Packet validity
    validity = rand() >= imu_params.packet_invalid_prob;

    % Build the output structs
    accelBusData.accel = accelerationSimulated;
    accelBusData.covariance = imu_params.acovariance;
    accelBusData.timestamp = timestamp;
    accelBusData.validity = validity;
    accelBus = accelBusData;

    angvelBusData.angvel = angularVelocitySimulated;
    angvelBusData.covariance = imu_params.gcovariance;
    angvelBusData.timestamp = timestamp;
    angvelBusData.validity = validity;
    angvelBus = angvelBusData;

end
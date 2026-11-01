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
    
    % Retrieve simulation parameters
    persistent imu_params;
    if isempty(imu_params)
        P = nav_sensor_params();   % defaults there; tests can override
        imu_params = P.imu;
    end
    
    timestamp = simtime;

    angularVelocityTrue = trueState(7:9);
    quaternionTrue = trueState(10:13);

    rotationInertialToBody = quaternion_to_rotation(quaternionTrue)';
    
    % Gravity in inertial NED coordinates
    gravityInertial = [0; 0; imu_params.gravity_scalar];
    
    % True specific force in body coordinates
    specificForceBodyTrue = rotationInertialToBody * (accelerationTrue(:) - gravityInertial);
    

    accelerometerNoise = imu_params.accel_noise_std .* randn(3, 1);
    gyroNoise = imu_params.gyro_noise_std .* randn(3,1);
    
    accelerationSimulated = specificForceBodyTrue + accelerometerBiasTrue(:) + accelerometerNoise;
    angularVelocitySimulated = angularVelocityTrue + gyroBiasTrue(:) + gyroNoise;
    angularVelocitySimulated = angularVelocitySimulated(1:3,1); % Simulink gets touchy about dims
    
    
    % Packet validity
    validity = rand() >= imu_params.packet_invalid_prob;

    % Build the output structs
    accelBus.accel = accelerationSimulated;
    accelBus.covariance = imu_params.acovariance;
    accelBus.timestamp = timestamp;
    accelBus.validity = validity;

    angvelBus.angvel = angularVelocitySimulated;
    angvelBus.covariance = imu_params.gcovariance;
    angvelBus.timestamp = timestamp;
    angvelBus.validity = validity;

end
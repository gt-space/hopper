function [accelerationSimulated, angularVelocitySimulated, timestamp, validity] = ...
    imu_model(trueState, accelerationTrue, accelerometerBiasTrue, gyroBiasTrue, simtime)
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
    %   accel      - Measured body-frame specific force (3x1)
    %   angvel     - Measured body-frame angular rate (3x1)
    %   timestamp  - Measurement timestamp
    %   validity   - Measurement validity flag
    
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

end
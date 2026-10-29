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
        imu_st = 1/500;
        imu_params.accel_noise_density = [880e-6; 880e-6; 732e-6];   % (m/s^2)/sqrt(Hz), x y z
        imu_params.gyro_noise_density = deg2rad([6.1e-3; 6.1e-3; 7.0e-3]); % (rad/s)/sqrt(Hz), x y z
        imu_params.accel_bias_stability = [125e-6; 125e-6; 134e-6];  % In-run bias stability (m/s^2)
        imu_params.gyro_bias_stability = deg2rad([7.5; 8.1; 4.9])/3600; % In-run bias stability (rad/s)
        imu_params.accel_noise_std = [880e-6; 880e-6; 732e-6] * sqrt(1/(2*imu_st)); % Accelerometer noise standard deviation (m/s^2)
        imu_params.acovariance = diag([880e-6; 880e-6; 732e-6].^2 * (1/(2*imu_st)));
        imu_params.gyro_noise_std = deg2rad([6.1e-3; 6.1e-3; 7.0e-3]) * sqrt(1/(2*imu_st)); % Gyroscope noise standard deviation (rad/s)
        imu_params.gcovariance = diag(deg2rad([6.1e-3; 6.1e-3; 7.0e-3].^2) * (1/(2*imu_st)));
        imu_params.gravity_scalar = 9.81;                % Gravity acceleration
        imu_params.packet_invalid_prob = 0.0;            % Set to 0 for now as requested
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
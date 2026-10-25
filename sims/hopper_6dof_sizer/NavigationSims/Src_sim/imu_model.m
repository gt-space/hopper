<<<<<<< HEAD:sims/hopper_6dof_sizer/NavigationSims/Src/imu_model.m
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
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imu_model.m
%
% Function to output the IMU measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
%
% Author: Aditiya Babbar
% Date: 10/11/26
%
% Inputs:
%   x: state vector (20x1)
%   a: acceleration vector (3x1)
%   accelBias: accelerometer measurement bias in body frame (m/s^2)
%   gyroBias: gyroscope measurement bias in body frame (rad/s)
%   simtime: time elapsed (s)
%
% Outputs:
%   accelBus: struct defined in the Navigation Dictionary
%   angvelBus: struct defined in the Navigation Dictionary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [accelBus, angvelBus] = imu_model(x, a, accelBias, gyroBias, simtime)

% Build empty output structs from the dictionary definition
accelBus = Simulink.Bus.createMATLABStruct('accelBus');
angvelBus = Simulink.Bus.createMATLABStruct('angvelBus');

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent imu_params;
if isempty(imu_params)
    try
        imu_params = evalin('base', 'imu_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        imu_params.accel_noise_std = [0; 0; 0];
        imu_params.gyro_noise_std = [0; 0; 0];
        imu_params.acovariance = eye(3)*[0.05; 0.05; 0.05].^2;
        imu_params.gcovariance = eye(3)*[0.005; 0.005; 0.005].^2;
        imu_params.gravity_scalar = 9.81; % Example body frame gravity component
        imu_params.packet_invalid_prob = 0.0;
>>>>>>> ef24f709a8dbddd6c79e22cd08112e04c7f92083:sims/hopper_6dof_sizer/NavigationSims/Src_sim/imu_model.m
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





<<<<<<< HEAD:sims/hopper_6dof_sizer/NavigationSims/Src/imu_model.m




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % imu_model.m
% %
% % Function to output the IMU measurement data with noise and bias applied.
% % Depends on filter_sim_initializer being ran first, which has bias/noise
% % params.
% %
% % Author: Aditiya Babbar
% % Date: 10/11/26
% %
% % Inputs:
% %   x: state vector (20x1)
% %   a: acceleration vector (3x1)
% %   accelBias: accelerometer measurement bias in body frame (m/s^2)
% %   gyroBias: gyroscope measurement bias in body frame (rad/s)
% %   simtime: time elapsed (s)
% %
% % Outputs:
% %   accel: acceleration vector (3x1) (m/s^2)
% %   angvel: angular velocity vector (3x1) (radians/s)
% %   timestamp: time elapsed (s)
% %   validity: boolean (true if valid measurement)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [accel, angvel, timestamp, validity] ...
%             = imu_model(x, a, accelBias, gyroBias, simtime)
% 
% % Retrieve parameters from the workspace (set by filter_sim_initializer)
% persistent imu_params;
% if isempty(imu_params)
%     try
%         imu_params = evalin('base', 'imu_params');
%     catch
%         % Fallback defaults if initializer hasn't populated workspace yet
%         imu_params.accel_noise_std = [0; 0; 0];
%         imu_params.gyro_noise_std = [0; 0; 0];
%         imu_params.gravity_scalar = 9.81; % Example body frame gravity component
%         imu_params.packet_invalid_prob = 0.0;
%     end
% end
% 
% % 1. Forward Timestamp
% timestamp = simtime;
% 
% % 2. Accelerometer Measurement: A(q)*(a_true + ba + na) - g_vec
% % Adjust gravity subtraction based on your coordinate convention/body frame definition
% g_vec = [0; 0; -g]; % accel reading in NED frame when IMU is still relative to ground
% 
% accel_ideal = a - g_radial;
% accel_noise = imu_params.accel_noise_std .* randn(3, 1);
% accel = accel_ideal + accelBias + accel_noise;
% 
% % 3. Gyroscope Measurement: w_true + bg + ng
% % NOTE: Update indices (e.g., x(10:12)) to match your specific 20x1 state vector layout for angular velocity
% omega_true = x(7:9);
% gyro_noise = imu_params.gyro_noise_std .* randn(3, 1);
% angvel = omega_true + gyroBias + gyro_noise;
% 
% % 4. Validity Flag Assignment
% if rand() < imu_params.packet_invalid_prob
%     validity = false;
% else
%     validity = true;
% end
% end
=======
% 4. Validity Flag Assignment
if rand() < imu_params.packet_invalid_prob
    validity = false;
else
    validity = true;
end

% Build the output structs
accelBus.accel = accel;
accelBus.covariance = imu_params.acovariance;
accelBus.timestamp = timestamp;
accelBus.validity = validity;

angvelBus.angvel = anvel;
angvelBus.covariance = imu_params.gcovariance;
angvelBus.timestamp = timestamp;
angvelBus.validity = validity;

end
>>>>>>> ef24f709a8dbddd6c79e22cd08112e04c7f92083:sims/hopper_6dof_sizer/NavigationSims/Src_sim/imu_model.m

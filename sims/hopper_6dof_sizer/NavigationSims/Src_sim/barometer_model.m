function [height, timestamp, validity] = barometer_model(trueState, baroBiasTrue, simtime)
    % Barometer measurement simulation.
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

    persistent baro_params;

    if isempty(baro_params)
        try
            baro_params = evalin('base', 'baro_params');
        catch
            baro_params.noise_std = 0.1;
            baro_params.packet_invalid_prob = 0.0;
        end
    end

    % validate inputs
    assert(numel(trueState) == 20);
    assert(isscalar(baroBiasTrue));
    assert(isscalar(simtime));
    
    timestamp = simtime;
    
    % NED position -> altitude relative to NED origin
    heightTrue = -trueState(3);
    
    % Measurement bias and Gaussian white noise
    baroNoise = baro_params.noise_std * randn();
    
    height = heightTrue + baroBiasTrue + baroNoise;
    


    validity = rand() >= baro_params.packet_invalid_prob;

end
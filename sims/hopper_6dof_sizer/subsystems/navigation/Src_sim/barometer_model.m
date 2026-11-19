% Barometer measurement simulation.
%
% Inputs:
%   trueState  - True vehicle state (20x1)
%   baroBiasTrue - True barometer height bias (m)
%   simtime    - Simulation time (s)
%   measurementNoise - Optional additive Sensors measurement noise (26x1)
%
% Outputs:
%   baroOut    - struct defined in the Navigation Dictionary
function baroOut = barometer_model(trueState, baroBiasTrue, simtime, measurementNoise)
    
    % Simulink gets touchy about dims
    trueState = trueState(1:20,1);

    % Load the baro_params struct
    persistent baro_params L;
    if isempty(baro_params)
        P = nav_sensor_params();   % defaults there; tests can override
        baro_params = P.baro;
        L = sensor_layout();
    end
    
    % Forward timestamp for now
    timestamp = simtime;
    
    % NED position -> altitude relative to NED origin
    heightTrue = -trueState(3,1);
    
    % Measurement bias and Gaussian white noise
    hasInjectedNoise = nargin >= 4 && ~isempty(measurementNoise);
    if hasInjectedNoise
        baroNoise = measurementNoise(L.BARO);
    else
        baroNoise = baro_params.noise_std * randn();
    end
    
    % Compute the measured height and assign validity
    height = heightTrue + baroBiasTrue + baroNoise;
    if hasInjectedNoise
        validity = true;
    else
        validity = rand() >= baro_params.packet_invalid_prob;
    end

    % Build the output struct
    baroOut.height = height;
    baroOut.covariance = baro_params.covariance;
    baroOut.timestamp = timestamp;
    baroOut.validity = validity;

end
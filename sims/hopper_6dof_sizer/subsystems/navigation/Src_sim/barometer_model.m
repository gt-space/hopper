% Barometer measurement simulation.
%
% Inputs:
%   trueState  - True vehicle state (20x1)
%   baroBiasTrue - True barometer height bias (m)
%   simtime    - Simulation time (s)
%
% Outputs:
%   baroOut    - struct defined in the Navigation Dictionary
function baroOut = barometer_model(trueState, baroBiasTrue, simtime)
    
    % Simulink gets touchy about dims
    trueState = trueState(1:20,1);

    % Load the baro_params struct
    persistent baro_params;
    if isempty(baro_params)
        P = nav_sensor_params();   % defaults there; tests can override
        baro_params = P.baro;
    end
    
    % Forward timestamp for now
    timestamp = simtime;
    
    % NED position -> altitude relative to NED origin
    heightTrue = -trueState(3,1);
    
    % Measurement bias and Gaussian white noise
    baroNoise = baro_params.noise_std * randn();
    
    % Compute the measured height and assign validity
    height = heightTrue + baroBiasTrue + baroNoise;
    validity = rand() >= baro_params.packet_invalid_prob;

    % Build the output struct
    baroOut.height = height;
    baroOut.covariance = baro_params.covariance;
    baroOut.timestamp = timestamp;
    baroOut.validity = validity;

end
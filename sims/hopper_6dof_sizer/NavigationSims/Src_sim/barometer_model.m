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
    
    % Load the baro_params struct
    persistent baro_params;
    if isempty(baro_params)
        baro_params.noise_std = 0.1;            % Noise standard deviation (m)
        baro_params.covariance = 0.1^2;         % Filter measurement covariance
        baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now
        baro_params.bias = 0.5;                 % Placeholder
    end

    % validate inputs
    assert(numel(trueState) == 20);
    assert(isscalar(simtime));
    
    timestamp = simtime;
    
    % NED position -> altitude relative to NED origin
    heightTrue = -trueState(3);
    
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
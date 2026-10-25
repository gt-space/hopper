% Barometer measurement simulation.
%
% Inputs:
%   trueState  - True vehicle state (20x1)
%   baroBiasTrue - True barometer height bias (m)
%   simtime    - Simulation time (s)
%
% Outputs:
%   baroBus    - struct defined in the Navigation Dictionary
function baroBus = barometer_model(trueState, baroBiasTrue, simtime)

    % Build an empty output struct from the dictionary definition
    persistent baroBusData;
    if isempty(baroBusData)
        dictionary = Simulink.data.dictionary.open('NavigationDictionary.sldd');
        baroBusData = Simulink.Bus.createMATLABStruct('baroBus',[],1,dictionary);
    end
    
    % Load the baro_params struct
    persistent baro_params;
    if isempty(baro_params)
        try
            baro_params = evalin('base', 'baro_params');
        catch
            baro_params.noise_std = 0.1;
            baro_params.covariance = 0.1^2;         % Filter measurement covariance
            baro_params.bias = 0.5; % placeholder
            baro_params.packet_invalid_prob = 0.0;
        end
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
    baroBusData.height = height;
    baroBusData.covariance = baro_params.covariance;
    baroBusData.timestamp = timestamp;
    baroBusData.validity = validity;
    baroBus = baroBusData;

end
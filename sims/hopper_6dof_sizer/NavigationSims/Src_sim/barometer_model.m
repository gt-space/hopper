<<<<<<< HEAD:sims/hopper_6dof_sizer/NavigationSims/Src/barometer_model.m
function [height, timestamp, validity] = barometer_model(trueState, baroBiasTrue, simtime)

    persistent baro_params;

    if isempty(baro_params)
        try
            baro_params = evalin('base', 'baro_params');
        catch
            baro_params.noise_std = 0.1;
            baro_params.packet_invalid_prob = 0.0;
        end
=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% barometer_model.m
%
% Function to output the baro measurement data with noise and bias applied.
% Depends on filter_sim_initializer being ran first, which has bias/noise
% params.
%
% Author: Aditiya Babbar
% Date: 10/11/26
%
% Inputs:
%   x: state vector (20x1)
%   simtime: time elapsed (s)
%
% Outputs:
%   baroBus: struct defined in the Navigation Dictionary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function baroBus = barometer_model(x, time)

% Build an empty output struct from the dictionary definition
baroBus = Simulink.Bus.createMATLABStruct('baroBus');

% Retrieve parameters from the workspace (set by filter_sim_initializer)
persistent baro_params;
if isempty(baro_params)
    try
        baro_params = evalin('base', 'baro_params');
    catch
        % Fallback defaults if initializer hasn't populated workspace yet
        baro_params.noise_std = 0.1;            % Height noise std deviation (m)
        baro_params.covariance = 0.1^2;
        baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now
        baro_params.bias = 0;
>>>>>>> ef24f709a8dbddd6c79e22cd08112e04c7f92083:sims/hopper_6dof_sizer/NavigationSims/Src_sim/barometer_model.m
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
    


<<<<<<< HEAD:sims/hopper_6dof_sizer/NavigationSims/Src/barometer_model.m
    validity = rand() >= baro_params.packet_invalid_prob;

end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % barometer_model.m
% %
% % Function to output the baro measurement data with noise and bias applied.
% % Depends on filter_sim_initializer being ran first, which has bias/noise
% % params.
% %
% % Author: Aditiya Babbar
% % Date: 10/11/26
% %
% % Inputs:
% %   x: state vector (20x1)
% %   baroBias: barometer measurement bias (m)
% %   simtime: time elapsed (s)
% %
% % Outputs:
% %   height: height from ground (m)
% %   timestamp: time elapsed (s)
% %   validity: boolean (true if valid measurement)
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [height, timestamp, validity] = barometer_model(x, baroBias, time)
% 
% % Retrieve parameters from the workspace (set by filter_sim_initializer)
% persistent baro_params;
% if isempty(baro_params)
%     try
%         baro_params = evalin('base', 'baro_params');
%     catch
%         % Fallback defaults if initializer hasn't populated workspace yet
%         baro_params.noise_std = 0.1;            % Height noise std deviation (m)
%         baro_params.packet_invalid_prob = 0.0;  % Set to 0 for now
%     end
% end
% 
% % 1. Forward Timestamp
% timestamp = time;
% 
% % 2. Compute Height Measurement with Bias and Noise
% 
% % Adjust index x(3) to match whatever index represents your true height/altitude in the 20x1 state vector
% true_height = x(3);
% baro_noise = baro_params.noise_std * randn();
% height = true_height + baroBias + baro_noise;
% 
% % 3. Validity Flag Assignment
% if rand() < baro_params.packet_invalid_prob
%     validity = false;
% else
%     validity = true;
% end
% end
=======
% Adjust index x(3) to match whatever index represents your true height/altitude in the 20x1 state vector
true_height = x(3);
baro_noise = baro_params.noise_std * randn();
height = true_height + baro_params.bias + baro_noise;

% 3. Validity Flag Assignment
validity = rand() > baro_params.packet_invalid_prob;

% Build the output struct
baroBus.height = height;
baroBus.covariance = baro_params.covariance;
baroBus.timestamp = timestamp;
baroBus.validity = validity;

end
>>>>>>> ef24f709a8dbddd6c79e22cd08112e04c7f92083:sims/hopper_6dof_sizer/NavigationSims/Src_sim/barometer_model.m

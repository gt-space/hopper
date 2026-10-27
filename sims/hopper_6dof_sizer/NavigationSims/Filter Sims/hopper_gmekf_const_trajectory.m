% hopper_gmekf_const_trajectory.m
% 
% Simulation of a constant trajectory to verify filter converges to the
% steady vehicle state.
% 
% Author: David Meisinger
% Date: 10/25/26

% First run the initializer to load noise, frame params, etc.
% Note that tuning Qc lives inside the initializer.
run('./filter_sim_initializer.m');

% Set up time span
tf = 100;       % stop time in seconds
Ts = 1/500;     % time increment in seconds
timespan = (0:Ts:tf)';
N = numel(timespan); % number of sample points

% Make sure the model is on the MATLAB path
load_system("Hopper_GMEKF_Simulation");

%% Configure simulation timing
set_param("Hopper_GMEKF_Simulation", ...
    'StartTime',       num2str(0), ...
    'StopTime',        num2str(tf), ...
    'SolverType',      'Fixed-step', ...
    'FixedStep',       num2str(Ts));
set_param("Hopper_GMEKF_Simulation", 'Solver', 'FixedStepDiscrete');

%% Define filter initialization (override filter_initialization if desired)
% hovering at 1m above ground with zero attitude and no biases.
x0 = [0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]';
% P0 = ...
% Qc = ...
% massProps0 = ...
gravity = 9.81;

%% Generate input time signals
accelTruth = zeros(N, 3);   % [ax ay az]
torqueTruth = zeros(N, 3);  % [Tx Ty Tz]
stateTruth = zeros(N, 20);  % 20-state vector

% Convert to Simulink time-series inputs
accelTruth_ts = timeseries(accelTruth, timespan);
torqueTruth_ts = timeseries(torqueTruth, timespan);
stateTruth_ts = timeseries(stateTruth, timespan);

%% Run simulation
simIn = Simulink.SimulationInput("Hopper_GMEKF_Simulation");

% Pass constant workspace variables explicitly to the simulation
simIn = simIn.setVariable('x0', x0);
simIn = simIn.setVariable('P0', P0);
simIn = simIn.setVariable('timestep', Ts);
simIn = simIn.setVariable('gravity', gravity);

out = sim(simIn);

%% Make plots
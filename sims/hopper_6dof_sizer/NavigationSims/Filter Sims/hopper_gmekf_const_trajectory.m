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
load_system("Hopper_GMEKF");
load_system("Hopper_GMEKF_Simulation");

% Open the data dictionary connection
sldd_name = 'NavigationDictionary.sldd';
dictObj = Simulink.data.dictionary.open(sldd_name);

% Access the Design Data section where ValueTypes live
designData = getSection(dictObj, 'Design Data');

%% Configure simulation timing
set_param("Hopper_GMEKF_Simulation", ...
    'StartTime',       num2str(0), ...
    'StopTime',        num2str(tf), ...
    'SolverType',      'Fixed-step', ...
    'FixedStep',       num2str(Ts));
set_param("Hopper_GMEKF_Simulation", 'Solver', 'FixedStepDiscrete');
set_param('Hopper_GMEKF', 'SolverType', 'Fixed-step')

%% Define filter initialization (override filter_initialization if desired)
% hovering at 1m above ground with zero attitude and no biases.
x0 = [0 0 -0.9 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]';
% P0 = ...
% Qc = ...
% massProps0 = ...
g = 9.81;

%% Generate input time signals
mass = 130;
inertia_tensor = [62.28 0 0;
                  0 62.28 0;
                  0 0 1.048];
torques = [0; 0; 0];
acceleration_truth = [0; 0; 0];
accel_bias = [0.01; -0.01; 0];
gyro_bias = [-0.02; 0; 0.01];
baro_bias = 0;
state_truth = [0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]';

% Empty variables for populating
inertia_tensor_data = zeros(3,3,N);
torque_data = zeros(3,N);
acceleration_data = zeros(3,N);
accel_bias_data = zeros(3,N);
gyro_bias_data = zeros(3,N);
baro_bias = zeros(1,N);
state_truth_data = zeros(20,N);

% Populate vector and matrix variables
for idx = 1:N
    inertia_tensor_data(:,:,idx) = inertia_tensor;
    torque_data(:,idx) = torques;
    acceleration_data(:,idx) = acceleration_truth;
    accel_bias_data(:,idx) = accel_bias;
    gyro_bias_data(:,idx) = gyro_bias;
    state_truth_data(:,idx) = state_truth;
end

% Convert to Simulink time-series inputs
stateTruth_ts = timeseries(state_truth_data, timespan);
mass_ts = timeseries(mass, timespan);
inertiaTensor_ts = timeseries(inertia_tensor_data, timespan);
torque_ts = timeseries(torque_data, timespan);
accelerationTruth_ts = timeseries(acceleration_data, timespan);
accelBias_ts = timeseries(accel_bias_data, timespan);
gyroBias_ts = timeseries(gyro_bias_data, timespan);
baroBias_ts = timeseries(baro_bias, timespan);

%% Run simulation
simIn = Simulink.SimulationInput("Hopper_GMEKF_Simulation");

% Pass constant workspace variables explicitly to the simulation
simIn = simIn.setVariable('state_truth', stateTruth_ts);
simIn = simIn.setVariable('acceleration_truth', accelerationTruth_ts);
simIn = simIn.setVariable('gravity_in', g);
simIn = simIn.setVariable('x0', x0);
simIn = simIn.setVariable('P0', P0);
simIn = simIn.setVariable('timestep_in', Ts);
simIn = simIn.setVariable('torque_in', torque_ts);
simIn = simIn.setVariable('inertia_tensor', inertiaTensor_ts);
simIn = simIn.setVariable('mass', mass_ts);
simIn = simIn.setVariable('accel_bias', accelBias_ts);
simIn = simIn.setVariable('gyro_bias', gyroBias_ts);
simIn = simIn.setVariable('baro_bias', baroBias_ts);

% Run
out = sim(simIn);

%% Make plots
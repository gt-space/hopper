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
tf = 10;        % stop time in seconds
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
x0 = [0 0 -0.5 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]';
g = 9.81;

%% Generate input time signals
mass = 130*ones(1,N);
inertia_tensor = [62.28 0 0;
                  0 62.28 0;
                  0 0 1.048];
torques = [0; 0; 0];
acceleration_truth = [0; 0; 0];
accel_bias = [0.01; -0.01; 0];
gyro_bias = [-0.02; 0; 0.01];
baro_bias = 0*ones(1,N);
state_truth = [0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]';

% Empty variables for populating
inertia_tensor_data = zeros(3,3,N);
torque_data = zeros(3,N);
acceleration_data = zeros(3,N);
accel_bias_data = zeros(3,N);
gyro_bias_data = zeros(3,N);
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
model = "Hopper_GMEKF_Simulation";

% Get the root Inports in port order
inports = find_system(model,'SearchDepth',1,'BlockType','Inport');
portNum = cellfun(@(b) str2double(get_param(b,'Port')), inports);
[~,order] = sort(portNum);
inports = inports(order);
inNames = get_param(inports,'Name');

% Build a timeseries for every Inport, keyed by Inport name
sig = struct();
sig.state_truth        = stateTruth_ts;        % 20xN
sig.acceleration_truth = accelerationTruth_ts; % 3xN
sig.gravity_in         = timeseries(g*ones(1,N), timespan);
sig.x0                 = timeseries(repmat(x0, 1, N), timespan);   % 20xN
sig.P0                 = timeseries(repmat(P0, 1, 1, N), timespan);
sig.timestep_in        = timeseries(Ts*ones(N,1), timespan);
sig.torque_in          = torque_ts;            % 3xN
sig.inertia_tensor     = inertiaTensor_ts;     % 3x3xN
sig.mass               = mass_ts;
sig.accel_bias         = accelBias_ts;
sig.gyro_bias          = gyroBias_ts;
sig.baro_bias          = baroBias_ts;

% 3. Assemble the Dataset in port order
ds = Simulink.SimulationData.Dataset;
for k = 1:numel(inNames)
    ds = ds.addElement(sig.(inNames{k}), inNames{k});
end

simIn = Simulink.SimulationInput(model);
simIn = simIn.setExternalInput(ds);
% simIn = simIn.setModelParameter('LoadExternalInput','on');
% simIn = simIn.setVariable('Qc', qc);   % only things blocks resolve by name stay as variables
out = sim(simIn);

%% Make plots
xout = out.yout{1}.Values.Data;
Pout = out.yout{2}.Values.Data;
tout = out.tout;

posx = squeeze(xout(1,1,:));
posy = squeeze(xout(2,1,:));
posz = squeeze(xout(3,1,:));
velx = squeeze(xout(4,1,:));
vely = squeeze(xout(5,1,:));
velz = squeeze(xout(6,1,:));
omegax = squeeze(xout(7,1,:));
omegay = squeeze(xout(8,1,:));
omegaz = squeeze(xout(9,1,:));
q0 = squeeze(xout(10,1,:));
q1 = squeeze(xout(11,1,:));
q2 = squeeze(xout(12,1,:));
q3 = squeeze(xout(13,1,:));

figure(1)
plot(tout, posx, tout, posy, tout, posz);
xlabel('Time (s)')
ylabel('Position (m)');
legend('X', 'Y', 'Z');
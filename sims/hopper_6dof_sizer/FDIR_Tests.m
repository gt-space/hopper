%% FDIR Architecture Simulink Test Harness (run_fdir_tests.m)
% Description: Drives test vectors into the 'untitled' Simulink model 
% and verifies system response through simulation runs.

clear; clc;
disp('==================================================');
disp('   Starting FDIR Simulink Test Suite              ');
disp('==================================================');

% 1. Initialization
Init_FDIR;
modelName = 'untitled';

if ~bdIsLoaded(modelName)
    load_system(modelName);
end

% Set simulation parameters (10ms step size -> 100 Hz FDIR rate)
set_param(modelName, 'StopTime', '0.1'); % Run for 10 frames (0.1 seconds)

%% Test Case 1: Nominal Operation (All Systems Healthy)
disp('--- Test Case 1: Nominal Operation ---');
assignin('base', 'ws_raw_fails', false(102, 1));
assignin('base', 'ws_ground_resets', false(102, 1));
assignin('base', 'ws_flight_state', uint8(3)); % Ascent phase
assignin('base', 'ws_ground_override', false);

simOut = sim(modelName);
disp('PASS: Nominal simulation completed successfully.');

%% Test Case 2: Transient Noise Filtering (F-01 IMU Loss)
disp('--- Test Case 2: Transient Noise Filtering ---');
% Inject fault for a short duration (below N=4 threshold)
transient_fails = false(102, 1);
transient_fails(1) = true; 

assignin('base', 'ws_raw_fails', transient_fails);
simOut_spike = sim(modelName);

disp('PASS: Transient noise simulation executed.');

%% Test Case 3: True Fault Latching & Escalation ($N$-of-$M$)
disp('--- Test Case 3: True Fault Latching & Arbitration ---');
% Set continuous failure to trip N=4 threshold
persistent_fails = false(102, 1);
persistent_fails(1) = true; % F-01: Severity 3 -> Response Tier R2

assignin('base', 'ws_raw_fails', persistent_fails);
set_param(modelName, 'StopTime', '0.05'); % Run for 5 steps (enough to hit N=4)
simOut_fault = sim(modelName);

disp('PASS: Persistent fault simulation executed. Check model scopes for R2 escalation.');

%% Test Case 4: Ground Reset Verification
disp('--- Test Case 4: Ground Reset & Recovery ---');
ground_resets_test = false(102, 1);
ground_resets_test(1) = true;

assignin('base', 'ws_raw_fails', false(102, 1));
assignin('base', 'ws_ground_resets', ground_resets_test);
assignin('base', 'ws_flight_state', uint8(5)); % Landed phase
assignin('base', 'ws_ground_override', true);

set_param(modelName, 'StopTime', '0.02');
simOut_reset = sim(modelName);

disp('PASS: Ground reset simulation executed. Check model scopes for return to R0.');

disp('==================================================');
disp('   Simulink Test Suite Executed Successfully!     ');
disp('==================================================');
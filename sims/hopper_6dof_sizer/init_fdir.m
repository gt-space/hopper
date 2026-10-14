%% FDIR Executive Initialization Script (init_fdir.m)
% Description: Initializes workspace variables, configuration constants,
% enumerated types, and Simulink Bus objects for the lander FDIR architecture.

clear; clc;
disp('Initializing FDIR Architecture for Landed/Flight Operations...');

%% 1. Configuration Parameters
FDIR_CONFIG.M_WINDOW_SIZE = 10;      % Rolling window size (M frames)
FDIR_CONFIG.N_FAIL_THRESH = 4;       % Failure threshold (N out of M)
FDIR_CONFIG.SAMPLE_RATE_SEC = 0.01;  % Execution rate: 10 ms (100 Hz for 100ms latency proof)

%% 2. Flight Phase Enumeration Check
% Ensure FlightPhase.m class file is in the current working folder
if exist('FlightPhase', 'class') ~= 8
    warning('FlightPhase enum class not found in path. Ensure FlightPhase.m is saved in your working directory.');
else
    disp('-> FlightPhase enumeration loaded successfully.');
end

%% 3. Fault State Bus Definition
% Clear existing bus object from base workspace if present
if evalin('base', 'exist(''Fault_State_Bus'', ''var'')')
    evalin('base', 'clear Fault_State_Bus;');
end

elems(1) = Simulink.BusElement;
elems(1).Name = 'Fault_ID';
elems(1).DataType = 'uint8';
elems(1).Dimensions = 1;
elems(1).Description = 'Unique identifier for the sensor/monitor';

elems(2) = Simulink.BusElement;
elems(2).Name = 'Monitor_Valid';
elems(2).DataType = 'boolean';
elems(2).Dimensions = 1;
elems(2).Description = 'True if sensor data is fresh and valid';

elems(3) = Simulink.BusElement;
elems(3).Name = 'Raw_Fail';
elems(3).DataType = 'boolean';
elems(3).Dimensions = 1;
elems(3).Description = 'Immediate pass/fail evaluation result';

elems(4) = Simulink.BusElement;
elems(4).Name = 'Persistence_State';
elems(4).DataType = 'uint8';
elems(4).Dimensions = 1;
elems(4).Description = 'State machine status: 0=INACTIVE, 1=PENDING, 2=LATCHED';

elems(5) = Simulink.BusElement;
elems(5).Name = 'Fail_Count';
elems(5).DataType = 'uint8';
elems(5).Dimensions = 1;
elems(5).Description = 'Number of failures in current M-window';

elems(6) = Simulink.BusElement;
elems(6).Name = 'First_Fault_Timestamp';
elems(6).DataType = 'double';
elems(6).Dimensions = 1;
elems(6).Description = 'Simulation time when fault first latched';

% Register the Fault State Bus object in the MATLAB base workspace
Fault_State_Bus = Simulink.Bus;
Fault_State_Bus.Elements = elems;
assignin('base', 'Fault_State_Bus', Fault_State_Bus);
disp('-> Fault_State_Bus object created and assigned to workspace.');

%% 4. Response Feedback Bus Definition (Mode Manager Handshake)
if evalin('base', 'exist(''Response_Feedback_Bus'', ''var'')')
    evalin('base', 'clear Response_Feedback_Bus;');
end

fb_elems(1) = Simulink.BusElement;
fb_elems(1).Name = 'Response_Ack';
fb_elems(1).DataType = 'boolean';
fb_elems(1).Dimensions = 1;
fb_elems(1).Description = 'Set true when mode manager acknowledges response';

fb_elems(2) = Simulink.BusElement;
fb_elems(2).Name = 'Execution_Status';
fb_elems(2).DataType = 'uint8';
fb_elems(2).Dimensions = 1;
fb_elems(2).Description = '0: Pending, 1: Active, 2: Completed, 3: Failed';

% Register the Response Feedback Bus object in the base workspace
Response_Feedback_Bus = Simulink.Bus;
Response_Feedback_Bus.Elements = fb_elems;
assignin('base', 'Response_Feedback_Bus', Response_Feedback_Bus);
disp('-> Response_Feedback_Bus object created and assigned to workspace.');

disp('FDIR Initialization Complete and Ready for Simulink Model Execution.');

%% FDIR Command Bus Definition
if evalin('base', 'exist(''FDIR_Command_Bus'', ''var'')')
    evalin('base', 'clear FDIR_Command_Bus;');
end

cmd_elems(1) = Simulink.BusElement;
cmd_elems(1).Name = 'Mask_Target';
cmd_elems(1).DataType = 'uint8';

cmd_elems(2) = Simulink.BusElement;
cmd_elems(2).Name = 'Degrade_Multiplier';
cmd_elems(2).DataType = 'double';

cmd_elems(3) = Simulink.BusElement;
cmd_elems(3).Name = 'Abort_Land_Active';
cmd_elems(3).DataType = 'boolean';

cmd_elems(4) = Simulink.BusElement;
cmd_elems(4).Name = 'Nav_Fallback_Active';
cmd_elems(4).DataType = 'boolean';

cmd_elems(5) = Simulink.BusElement;
cmd_elems(5).Name = 'Cutoff_Active';
cmd_elems(5).DataType = 'boolean';

cmd_elems(6) = Simulink.BusElement;
cmd_elems(6).Name = 'Entry_Trigger';
cmd_elems(6).DataType = 'boolean';

FDIR_Command_Bus = Simulink.Bus;
FDIR_Command_Bus.Elements = cmd_elems;
assignin('base', 'FDIR_Command_Bus', FDIR_Command_Bus);
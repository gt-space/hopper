%% FDIR Executive Initialization Script (init_fdir.m)
% Description: Initializes workspace variables, configuration constants,
% enumerated types, and Simulink Bus objects for the lander FDIR architecture.
% Fully traceable to CoDR GNC fault specifications[cite: 2].

clear; clc;
disp('Initializing FDIR Architecture for Landed/Flight Operations...');

%% 1. Configuration Parameters
FDIR_CONFIG.M_WINDOW_SIZE = 10;      % Rolling window size (M frames)[cite: 2]
FDIR_CONFIG.N_FAIL_THRESH = 4;       % Failure threshold (N out of M)[cite: 2]
FDIR_CONFIG.SAMPLE_RATE_SEC = 0.01;  % Execution rate: 10 ms (100 Hz for 100ms latency proof)[cite: 2]

%% 2. Flight Phase Enumeration Check
if exist('FlightPhase', 'class') ~= 8
    warning('FlightPhase enum class not found in path. Ensure FlightPhase.m is saved in your working directory.');
else
    disp('-> FlightPhase enumeration loaded successfully.');
end

%% 3. Fault State Bus Definition
if evalin('base', 'exist(''Fault_State_Bus'', ''var'')')
    evalin('base', 'clear Fault_State_Bus;');
end

elems(1) = Simulink.BusElement;
elems(1).Name = 'Fault_ID';
elems(1).DataType = 'uint8';
elems(1).Dimensions = 1;

elems(2) = Simulink.BusElement;
elems(2).Name = 'Monitor_Valid';
elems(2).DataType = 'boolean';
elems(2).Dimensions = 1;

elems(3) = Simulink.BusElement;
elems(3).Name = 'Raw_Fail';
elems(3).DataType = 'boolean';
elems(3).Dimensions = 1;

elems(4) = Simulink.BusElement;
elems(4).Name = 'Persistence_State';
elems(4).DataType = 'uint8';
elems(4).Dimensions = 1;

elems(5) = Simulink.BusElement;
elems(5).Name = 'Fail_Count';
elems(5).DataType = 'uint8';
elems(5).Dimensions = 1;

elems(6) = Simulink.BusElement;
elems(6).Name = 'First_Fault_Timestamp';
elems(6).DataType = 'double';
elems(6).Dimensions = 1;

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

fb_elems(2) = Simulink.BusElement;
fb_elems(2).Name = 'Execution_Status';
fb_elems(2).DataType = 'uint8';
fb_elems(2).Dimensions = 1;

Response_Feedback_Bus = Simulink.Bus;
Response_Feedback_Bus.Elements = fb_elems;
assignin('base', 'Response_Feedback_Bus', Response_Feedback_Bus);
disp('-> Response_Feedback_Bus object created and assigned to workspace.');

%% 5. FDIR Command Bus Definition
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
disp('-> FDIR_Command_Bus object created and assigned to workspace.');

%% 6. Fault Configuration Table (Table-Driven GNC FDIR)
% Uniform struct array for code generation compliance, traceable to CoDR[cite: 2]
clear fault_config;
for i = 1:5
    fault_config(i).ID            = uint8(i);
    fault_config(i).N_Threshold   = uint8(4);
    fault_config(i).M_Window      = uint8(10);
    
    % Map response tiers traceable to CoDR GNC fault documentation[cite: 2]
    if i <= 2
        fault_config(i).Response_Tier = uint8(1); % R1: Degrade / Mask[cite: 2]
    elseif i <= 4
        fault_config(i).Response_Tier = uint8(2); % R2: Abort to Land[cite: 2]
    elseif i == 5
        fault_config(i).Response_Tier = uint8(3); % R3: Nav-Degraded Abort[cite: 2]
    else
        fault_config(i).Response_Tier = uint8(4); % R4: Cutoff & Safing[cite: 2]
    end
end

assignin('base', 'fault_config', fault_config);
assignin('base', 'FDIR_CONFIG', FDIR_CONFIG);
disp('-> fault_config struct array created and assigned to workspace.');

disp('FDIR Initialization Complete and Ready for Simulink Model Execution.');
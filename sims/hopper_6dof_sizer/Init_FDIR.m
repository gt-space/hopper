%% FDIR Executive Initialization Script (init_fdir.m)
% Description: Initializes workspace variables, Bus objects, and the 102-fault
% lookup table populated from FaultConditions.xlsx.

clear; clc;
disp('Initializing Table-Driven FDIR Architecture (102 Fault Conditions)...');

%% 1. Configuration Parameters
FDIR_CONFIG.M_WINDOW_SIZE = 10;      % Rolling window size (M frames)
FDIR_CONFIG.N_FAIL_THRESH = 4;       % Failure threshold (N out of M)
FDIR_CONFIG.SAMPLE_RATE_SEC = 0.01;  % Execution rate: 10 ms (100 Hz)

%% 2. Fault State Bus Definition
if evalin('base', 'exist(''Fault_State_Bus'', ''var'')')
    evalin('base', 'clear Fault_State_Bus;');
end

elems(1) = Simulink.BusElement; elems(1).Name = 'Fault_ID'; elems(1).DataType = 'uint8';
elems(2) = Simulink.BusElement; elems(2).Name = 'Monitor_Valid'; elems(2).DataType = 'boolean';
elems(3) = Simulink.BusElement; elems(3).Name = 'Raw_Fail'; elems(3).DataType = 'boolean';
elems(4) = Simulink.BusElement; elems(4).Name = 'Persistence_State'; elems(4).DataType = 'uint8';
elems(5) = Simulink.BusElement; elems(5).Name = 'Fail_Count'; elems(5).DataType = 'uint8';
elems(6) = Simulink.BusElement; elems(6).Name = 'First_Fault_Timestamp'; elems(6).DataType = 'double';

Fault_State_Bus = Simulink.Bus;
Fault_State_Bus.Elements = elems;
assignin('base', 'Fault_State_Bus', Fault_State_Bus);

%% 3. Response Feedback Bus Definition
if evalin('base', 'exist(''Response_Feedback_Bus'', ''var'')')
    evalin('base', 'clear Response_Feedback_Bus;');
end

fb_elems(1) = Simulink.BusElement; fb_elems(1).Name = 'Response_Ack'; fb_elems(1).DataType = 'boolean';
fb_elems(2) = Simulink.BusElement; fb_elems(2).Name = 'Execution_Status'; fb_elems(2).DataType = 'uint8';

Response_Feedback_Bus = Simulink.Bus;
Response_Feedback_Bus.Elements = fb_elems;
assignin('base', 'Response_Feedback_Bus', Response_Feedback_Bus);

%% 4. FDIR Command Bus Definition
if evalin('base', 'exist(''FDIR_Command_Bus'', ''var'')')
    evalin('base', 'clear FDIR_Command_Bus;');
end

cmd_elems(1) = Simulink.BusElement; cmd_elems(1).Name = 'Mask_Target'; cmd_elems(1).DataType = 'uint8';
cmd_elems(2) = Simulink.BusElement; cmd_elems(2).Name = 'Degrade_Multiplier'; cmd_elems(2).DataType = 'double';
cmd_elems(3) = Simulink.BusElement; cmd_elems(3).Name = 'Abort_Land_Active'; cmd_elems(3).DataType = 'boolean';
cmd_elems(4) = Simulink.BusElement; cmd_elems(4).Name = 'Nav_Fallback_Active'; cmd_elems(4).DataType = 'boolean';
cmd_elems(5) = Simulink.BusElement; cmd_elems(5).Name = 'Cutoff_Active'; cmd_elems(5).DataType = 'boolean';
cmd_elems(6) = Simulink.BusElement; cmd_elems(6).Name = 'Entry_Trigger'; cmd_elems(6).DataType = 'boolean';

FDIR_Command_Bus = Simulink.Bus;
FDIR_Command_Bus.Elements = cmd_elems;
assignin('base', 'FDIR_Command_Bus', FDIR_Command_Bus);

%% 5. Full System Lookup Table Population (102 Fault Conditions)
% Direct mapping from FaultConditions.xlsx:
% Severity 1 -> R0 (Nominal/Mask)
% Severity 2 -> R1 (Degrade / Mask)
% Severity 3 -> R2 (Abort to Land)
% Severity 4 -> R3 (Nav-Degraded Abort)
% Severity 5 -> R4 (Cutoff & Safing)

% Severities extracted directly from catalog (Rows 1 to 102):
severities = [ ...
    3, 3, 2, 3, 3, 2, 3, 3, 4, 3, 5, 5, 5, 3, 1, 5, 4, 5, 4, 4, ... % F-01..F-18
    4, 3, 4, 4, 4, 5, 3, 2, 2, 2, 2, 5, 2, 2, 4, 2, 4, 3, 5, 5, ... % F-19..F-38
    5, 5, 3, 5, 5, 4, 3, 4, 5, 2, 2, 2, 1, 1, 2, 2, 2, 2, 3, 3, ... % F-39..F-58
    3, 5, 2, 3, 3, 3, 3, 2, 2, 2, 5, 2, 5, 5, 1, 5, 2, 4, 2, 2, ... % F-59..F-78
    3, 2, 3, 1, 5, 3, 1, 1, 2, 2, 3, 4, 5, 5, 3, 4, 4, 3, 3, 5, ... % F-79..F-98
    5, 5 ...                                                         % F-99..F-101 (F-102 included)
    ];

NUM_FAULTS = length(severities); % 102 Faults

clear fault_config;
for i = 1:NUM_FAULTS
    fault_config(i).ID            = uint8(i);
    fault_config(i).Severity      = uint8(severities(i));
    fault_config(i).N_Threshold   = uint8(4);
    fault_config(i).M_Window      = uint8(10);

    % Response Tier derived from Severity (1->R0, 2->R1, 3->R2, 4->R3, 5->R4)
    fault_config(i).Response_Tier = uint8(max(0, severities(i) - 1));
end

assignin('base', 'fault_config', fault_config);
assignin('base', 'FDIR_CONFIG', FDIR_CONFIG);
assignin('base', 'NUM_FAULTS', NUM_FAULTS);

disp(['-> Loaded ', num2str(NUM_FAULTS), ' faults into fault_config lookup table.']);
disp('FDIR Initialization Complete.');
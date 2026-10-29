%% FDIR Executive Initialization Script (Init_FDIR.m)
% Description: Initializes workspace variables, Bus objects, and the 102-fault
% lookup table populated from FaultConditions.xlsx.
%
% Tuning and the fault catalog live in fdir_config.m so the monitors and the
% tests share one source of truth. This script only publishes them and the
% Bus objects to the base workspace. It deliberately does not clear the
% caller's workspace.

disp('Initializing Table-Driven FDIR Architecture (102 Fault Conditions)...');

[FDIR_CONFIG, fault_config, fault_map, check_config] = fdir_config();
NUM_FAULTS = FDIR_CONFIG.NUM_FAULTS;

%% 1. Fault State Bus Definition
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

%% 2. Response Feedback Bus Definition
if evalin('base', 'exist(''Response_Feedback_Bus'', ''var'')')
    evalin('base', 'clear Response_Feedback_Bus;');
end

fb_elems(1) = Simulink.BusElement; fb_elems(1).Name = 'Response_Ack'; fb_elems(1).DataType = 'boolean';
fb_elems(2) = Simulink.BusElement; fb_elems(2).Name = 'Execution_Status'; fb_elems(2).DataType = 'uint8';

Response_Feedback_Bus = Simulink.Bus;
Response_Feedback_Bus.Elements = fb_elems;
assignin('base', 'Response_Feedback_Bus', Response_Feedback_Bus);

%% 3. FDIR Command Bus Definition
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

%% 4. Publish configuration and lookup table
assignin('base', 'fault_config', fault_config);
assignin('base', 'FDIR_CONFIG', FDIR_CONFIG);
assignin('base', 'fault_map', fault_map);
assignin('base', 'check_config', check_config);

% Readable catalog for the workspace: ID, subsystem, failure condition,
% severity, tier and what that tier commands. Not used by the model.
try
    assignin('base', 'fault_catalog', fdir_fault_catalog());
catch catalogErr
    warning('Init_FDIR:CatalogUnavailable', ...
        'Could not load fdir_fault_catalog.csv (%s).', catalogErr.message);
end
assignin('base', 'NUM_FAULTS', NUM_FAULTS);

disp(['-> Loaded ', num2str(NUM_FAULTS), ' faults into fault_config lookup table.']);
disp('FDIR Initialization Complete.');

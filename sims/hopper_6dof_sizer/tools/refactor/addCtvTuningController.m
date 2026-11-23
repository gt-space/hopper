function addCtvTuningController()
%ADDCTVTUNINGCONTROLLER  Build subsystems/gnc/LQR_CTV_Tuning.slx.
%   Makes the CONTROLLER=lqr_ctv_tuning choice (hopper_6dof/GNC/Controller/
%   Controller_LQR_CTV_Tuning) a Subsystem Reference. The file is the
%   lqr_sched choice's own scheduler (Clock, the reference time grid, the
%   reference state and the nominal command, all copied from
%   Controller_LQR_Sched) reading its own gain tables: K1flat_ctv and
%   K2grid_ctv, which src/gnc/ctvTuningGains.m designs for the CTV.
%
%   The scheduler's inputs (state, clock) and its command pass through
%   Zero-Order Holds whose sample times, ctrl_ts_in and ctrl_ts_out, come
%   from ctrlRateWorkspace: inherited (every solver step) by default, or the
%   flight computer's period and latency.
%
%   Harmless to rerun: it deletes and rebuilds the file. The harness block
%   already references the file by name.
%
%   The file must be a subsystem file (BlockDiagramType 'subsystem'). Saved
%   as a model it can't be referenced: Simulink rejects the block with
%   "Invalid block diagram type 'Model'".

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % sims/hopper_6dof_sizer
file = fullfile(root, 'subsystems', 'gnc', 'LQR_CTV_Tuning.slx');
name = 'LQR_CTV_Tuning';

load_system('hopper_6dof');
src = 'hopper_6dof/GNC/Controller/Controller_LQR_Sched';   % the working choice

if bdIsLoaded(name), close_system(name, 0); end
if isfile(file), delete(file); end
new_system(name, 'Subsystem');
open_system(name);

% the choice's blocks: x, u, Clock, the constants and the scheduler reference
for b = get_param(src, 'Blocks')'
    add_block([src '/' b{1}], [name '/' b{1}], 'MakeNameUnique', 'off');
end
% and its wiring, by port number
for l = find_system(src, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'Line')'
    sh = get_param(l, 'SrcPortHandle');
    dh = get_param(l, 'DstPortHandle');
    if sh <= 0 || all(dh <= 0), continue; end
    sb = get_param(sh, 'Parent');
    si = find(get_param(sb, 'PortHandles').Outport == sh, 1);
    for d = dh(:)'
        db = get_param(d, 'Parent');
        di = find(get_param(db, 'PortHandles').Inport == d, 1);
        add_line(name, sprintf('%s/%d', get_param(sb, 'Name'), si), ...
            sprintf('%s/%d', get_param(db, 'Name'), di), 'autorouting', 'smart');
    end
end

% the tables this variant designs for itself
tables = struct('K1flat', 'K1flat_ctv', 'K2grid', 'K2grid_ctv');
for b = get_param(name, 'Blocks')'
    p = [name '/' b{1}];
    try
        v = get_param(p, 'Value');
    catch
        continue
    end
    if isfield(tables, v)
        set_param(p, 'Value', tables.(v));
    end
end

% The flight computer's timing (ctrlRateWorkspace): the scheduler reads
% the state and the clock at each control period's start (ctrl_ts_in) and
% its command lands ctrl_latency into the period (ctrl_ts_out). Both are
% -1 (inherited) by default, which leaves the holds as pass-throughs and
% the controller running every solver step.
insertHold(name, 'x', 'Ctrl_State_Hold', 'ctrl_ts_in', 'src');
insertHold(name, 'Clock', 'Ctrl_Clock_Hold', 'ctrl_ts_in', 'src');
insertHold(name, 'u', 'Ctrl_Command_Hold', 'ctrl_ts_out', 'dst');

save_system(name, file);
close_system(name, 0);
load_system(name);
if get_param(name, 'BlockDiagramType') ~= "subsystem"
    error('addCtvTuningController:type', '%s is a %s, not a subsystem file.', ...
        file, get_param(name, 'BlockDiagramType'));
end
ips = find_system(name, 'SearchDepth', 1, 'BlockType', 'Inport');
ops = find_system(name, 'SearchDepth', 1, 'BlockType', 'Outport');
close_system(name, 0);
if isempty(ips) || isempty(ops)
    error('addCtvTuningController:ports', '%s has no inport/outport at its top level.', file);
end
fprintf('%s: %d in, %d out, reads K1flat_ctv and K2grid_ctv\n', file, numel(ips), numel(ops));
end

function insertHold(sys, block, holdName, tsVar, side)
% Puts a Zero-Order Hold (sample time tsVar) on every line leaving block
% (side 'src') or on the line into block (side 'dst').
hold = [sys '/' holdName];
pos = get_param([sys '/' block], 'Position');
if side == "src", dx = 60; else, dx = -60; end
add_block('simulink/Discrete/Zero-Order Hold', hold, 'SampleTime', tsVar, ...
    'Position', [pos(1) + dx, pos(2), pos(1) + dx + 30, pos(2) + 24]);
ph = get_param([sys '/' block], 'PortHandles');
if side == "src"
    l = get_param(ph.Outport(1), 'Line');
    dsts = get_param(l, 'DstPortHandle');
    delete_line(l);
    hp = get_param(hold, 'PortHandles');
    add_line(sys, ph.Outport(1), hp.Inport(1), 'autorouting', 'smart');
    for d = dsts(:)'
        add_line(sys, hp.Outport(1), d, 'autorouting', 'smart');
    end
else
    l = get_param(ph.Inport(1), 'Line');
    s = get_param(l, 'SrcPortHandle');
    delete_line(l);
    hp = get_param(hold, 'PortHandles');
    add_line(sys, s, hp.Inport(1), 'autorouting', 'smart');
    add_line(sys, hp.Outport(1), ph.Inport(1), 'autorouting', 'smart');
end
end

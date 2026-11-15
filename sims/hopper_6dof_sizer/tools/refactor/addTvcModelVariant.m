function addTvcModelVariant()
%ADDTVCMODELVARIANT  TVC_MODEL variant axis in Actuators_TVC_RCS.
%   subsystems/actuators/Actuators_TVC_RCS.slx: the per-axis gimbal chain
%   (Saturation -> Transport Delay -> Quantizer -> Rate Limiter ->
%   Transfer Fcn, and the "1" blocks for yaw) moves into the Variant
%   Subsystem TVC_Actuator:
%     Second_Order  (VARIANT_TVC_MODEL == 1) the same blocks, unchanged,
%                   so the default flight is bit for bit the old one
%     Electromech   (VARIANT_TVC_MODEL >= 2) MATLAB Function TVC_EM: the
%                   Erhardt+Leimer voltage-to-gimbal-angle model
%                   (src/tvc/tvcEmUpdate.m) at the 1 kHz servo rate,
%                   parameters from the workspace struct TVCEM (tvcWorkspace)
%   TVC_Actuator takes the pitch and yaw commands and the thrust (the
%   actuator load), and its outputs feed the bending sum and modal path
%   exactly where the Transfer Fcn outputs did. Rerunning is harmless.

f = 'Actuators_TVC_RCS';
load_system(f);
vss = [f '/TVC_Actuator'];
if getSimulinkBlockHandle(vss) > 0
    fprintf('%s: already there\n', vss);
    close_system(f, 0);
    return
end

chain = {'Saturation', 'Transport Delay', 'Quantizer', 'Rate Limiter', 'Transfer Fcn'};
yawOf = @(b) [b '1'];
p0 = get_param([f '/Saturation'], 'Position');
p1 = get_param([f '/Transfer Fcn1'], 'Position');

% ---- the variant subsystem and its ports ----
add_block('simulink/Ports & Subsystems/Variant Subsystem', vss, ...
    'Position', [p0(1), p0(2), p1(3), p0(4) + 80]);
% empty the library template (its example choices and ports)
for b = find_system(vss, 'SearchDepth', 1, 'Type', 'block')'
    if ~strcmp(b{1}, vss), delete_block(b{1}); end
end
delete_line(find_system(vss, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line'));
ports = {'cmd_pitch', 'cmd_yaw', 'thrust'};
outs = {'delta_pitch', 'delta_yaw'};
for k = 1:3
    add_block('built-in/Inport', [vss '/' ports{k}], 'Position', [20, 40*k, 50, 40*k + 14]);
end
for k = 1:2
    add_block('built-in/Outport', [vss '/' outs{k}], 'Position', [600, 40*k, 630, 40*k + 14]);
end

% ---- choice 1: the existing chain ----
c1 = [vss '/Second_Order'];
add_block('built-in/SubSystem', c1, 'Position', [200, 30, 400, 110]);
set_param(c1, 'VariantControl', 'VARIANT_TVC_MODEL == 1');
addPorts(c1, ports, outs);
for ax = 1:2
    prev = [ports{ax} '/1'];
    for k = 1:numel(chain)
        src = chain{k};  if ax == 2, src = yawOf(src); end
        dst = [c1 '/' src];
        add_block([f '/' src], dst, 'Position', [100 + 120*k, 60*ax, 180 + 120*k, 60*ax + 30]);
        add_line(c1, prev, [src '/1'], 'autorouting', 'on');
        prev = [src '/1'];
    end
    add_line(c1, prev, [outs{ax} '/1'], 'autorouting', 'on');
end
add_block('built-in/Terminator', [c1 '/Unused thrust'], 'Position', [140, 200, 160, 220]);
add_line(c1, 'thrust/1', 'Unused thrust/1', 'autorouting', 'on');

% ---- choice 2: electromechanical model ----
c2 = [vss '/Electromech'];
add_block('built-in/SubSystem', c2, 'Position', [200, 150, 400, 230]);
set_param(c2, 'VariantControl', 'VARIANT_TVC_MODEL >= 2');
addPorts(c2, ports, outs);
add_block('built-in/Mux', [c2 '/Mux'], 'Inputs', '2', 'Position', [150, 40, 155, 110]);
add_block('built-in/Demux', [c2 '/Demux'], 'Outputs', '2', 'Position', [470, 40, 475, 110]);
em = [c2 '/TVC_EM'];
add_block('simulink/User-Defined Functions/MATLAB Function', em, 'Position', [240, 50, 380, 130]);
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', em);
ch.Script = strjoin({
    'function delta = TVC_EM(cmd, thrust, TVCEM)'
    '% Erhardt+Leimer TVC actuators, voltage to gimbal angle (src/tvc/tvcEmUpdate.m),'
    '% one call per 1 kHz servo period; cmd and delta are [pitch; yaw] in rad.'
    '% (tools/refactor/addTvcModelVariant.m)'
    'persistent S'
    'if isempty(S), S = tvcEmInit(); end'
    '[S, delta] = tvcEmUpdate(S, cmd, thrust, TVCEM);'
    'end'}, newline);
d = ch.find('-isa', 'Stateflow.Data', 'Name', 'TVCEM');
d.Scope = 'Parameter';
d.Tunable = false;
o = ch.find('-isa', 'Stateflow.Data', 'Name', 'delta');
o.Props.Array.Size = '[2 1]';
ch.SampleTime = '1e-3';
% One servo period of input delay breaks the algebraic loop through the
% thrust and the vehicle state (a MATLAB Function block is direct
% feedthrough); tvcWorkspace takes the period out of TVCEM.lat_steps so
% the command latency stays IN.tvc.actuator_latency.
add_block('built-in/UnitDelay', [c2 '/Cmd Delay'], 'SampleTime', '1e-3', 'Position', [185, 60, 215, 90]);
add_block('built-in/UnitDelay', [c2 '/Thrust Delay'], 'SampleTime', '1e-3', 'Position', [185, 150, 215, 180]);
add_line(c2, 'cmd_pitch/1', 'Mux/1', 'autorouting', 'on');
add_line(c2, 'cmd_yaw/1', 'Mux/2', 'autorouting', 'on');
add_line(c2, 'Mux/1', 'Cmd Delay/1', 'autorouting', 'on');
add_line(c2, 'Cmd Delay/1', 'TVC_EM/1', 'autorouting', 'on');
add_line(c2, 'thrust/1', 'Thrust Delay/1', 'autorouting', 'on');
add_line(c2, 'Thrust Delay/1', 'TVC_EM/2', 'autorouting', 'on');
add_line(c2, 'TVC_EM/1', 'Demux/1', 'autorouting', 'on');
add_line(c2, 'Demux/1', 'delta_pitch/1', 'autorouting', 'on');
add_line(c2, 'Demux/2', 'delta_yaw/1', 'autorouting', 'on');

% ---- rewire the top level ----
for ax = 1:2
    for k = 1:numel(chain)
        b = chain{k};  if ax == 2, b = yawOf(b); end
        lh = get_param([f '/' b], 'LineHandles');
        for l = [lh.Inport, lh.Outport]
            if l > 0, delete_line(l); end
        end
        delete_block([f '/' b]);
    end
end
add_line(f, 'TVC Pitch/1', 'TVC_Actuator/1', 'autorouting', 'on');
add_line(f, 'TVC Yaw/1', 'TVC_Actuator/2', 'autorouting', 'on');
add_line(f, 'Thrust /1', 'TVC_Actuator/3', 'autorouting', 'on');
add_line(f, 'TVC_Actuator/1', 'Sum5/1', 'autorouting', 'on');
add_line(f, 'TVC_Actuator/1', 'Transfer Fcn2/1', 'autorouting', 'on');
add_line(f, 'TVC_Actuator/2', 'Sum8/1', 'autorouting', 'on');
add_line(f, 'TVC_Actuator/2', 'Transfer Fcn3/1', 'autorouting', 'on');

save_system(f);
fprintf('%s: Second_Order (VARIANT_TVC_MODEL == 1) and Electromech (>= 2)\n', vss);
close_system(f, 0);
end

function addPorts(sys, ins, outs)
for k = 1:numel(ins)
    add_block('built-in/Inport', [sys '/' ins{k}], 'Position', [20, 60*k, 50, 60*k + 14]);
end
for k = 1:numel(outs)
    add_block('built-in/Outport', [sys '/' outs{k}], 'Position', [800, 60*k, 830, 60*k + 14]);
end
end

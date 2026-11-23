function addFswControllerChoice()
%ADDFSWCONTROLLERCHOICE  CONTROLLER=fsw: flight software flies the harness.
%   Adds the CONTROLLER=fsw choice to models/hopper_6dof.slx, so the
%   FSW bridge (fsw_bridge/) can generate the plant straight from the
%   harness instead of from a copied model:
%   - GNC/Controller gets an input u_fsw and an output x_fsw, and a choice
%     Controller_FSW (VARIANT_CONTROLLER == fsw's index in variantAxes, 5
%     since lqr_ctv_tuning took 4) that passes u_fsw through as
%     the command u and hands the controller's state x out as x_fsw.
%   - Both ports run up through GNC to a root Inport u_fsw (4 x 1: thrust
%     N, TVC pitch rad, TVC yaw rad, RCS torque N m; the choice takes it
%     fsw_latency into each fsw_period and holds it)
%     and a root Outport x_fsw (the 13-state the controller flies on).
%   - GNC/Controller propagates its variant conditions, so with any other
%     controller the new ports drop out of the compiled model and the
%     generated code: the other flights are unchanged.
%   Rerunning is harmless.

m = 'hopper_6dof';
load_system(m);
ctl = [m '/GNC/Controller'];
choice = [ctl '/Controller_FSW'];
if getSimulinkBlockHandle(choice) > 0
    fprintf('%s: Controller_FSW already there\n', m);
    close_system(m, 0);
    return
end

% ---- the choice inside the Controller variant subsystem ----
add_block('simulink/Sources/In1', [ctl '/u_fsw'], 'Position', blockPos(ctl, 'x', [0 80]));
add_block('simulink/Sinks/Out1', [ctl '/x_fsw'], 'Position', blockPos(ctl, 'u', [0 80]));
add_block('simulink/Ports & Subsystems/Subsystem', choice, ...
    'Position', blockPos(ctl, 'Controller_LQR_Sched_Legacy', [0 120]));
delete_line(choice, 'In1/1', 'Out1/1');
delete_block([choice '/In1']);
delete_block([choice '/Out1']);
add_block('simulink/Sources/In1', [choice '/x'], 'Position', [40 40 70 54]);
add_block('simulink/Sources/In1', [choice '/u_fsw'], 'Position', [40 100 70 114]);
add_block('simulink/Sinks/Out1', [choice '/u'], 'Position', [240 100 270 114]);
add_block('simulink/Sinks/Out1', [choice '/x_fsw'], 'Position', [240 40 270 54]);
% The flight software reads the state at the start of each control period
% and its command lands fsw_latency later (one plant step): the hold
% samples u_fsw at that offset and keeps it for the rest of the period
add_block('simulink/Discrete/Zero-Order Hold', [choice '/FSW_Period_Hold'], ...
    'SampleTime', '[fsw_period, fsw_latency]', 'Position', [130 95 165 119]);
add_line(choice, 'u_fsw/1', 'FSW_Period_Hold/1');
add_line(choice, 'FSW_Period_Hold/1', 'u/1');
add_line(choice, 'x/1', 'x_fsw/1');
% the index comes from variantAxes, never a literal: two branches once both
% claimed 4 (fsw and lqr_ctv_tuning)
k = find(variantAxes().CONTROLLER == "fsw");
assert(isscalar(k), 'addFswControllerChoice:axis', 'variantAxes has no single fsw choice');
set_param(choice, 'VariantControl', sprintf('VARIANT_CONTROLLER == %d', k));
set_param(ctl, 'PropagateVariantConditions', 'on');

% ---- through GNC ----
gnc = [m '/GNC'];
add_block('simulink/Sources/In1', [gnc '/u_fsw'], 'Position', blockPos(gnc, 'Controller', [-200 60], [30 14]));
add_block('simulink/Sinks/Out1', [gnc '/x_fsw'], 'Position', blockPos(gnc, 'u_cmd', [0 120]));
add_line(gnc, 'u_fsw/1', 'Controller/2', 'autorouting', 'smart');
add_line(gnc, 'Controller/2', 'x_fsw/1', 'autorouting', 'smart');

% ---- root ports ----
inp = [m '/u_fsw'];
add_block('simulink/Sources/In1', inp, 'Position', blockPos(m, 'GNC', [-220 -60], [30 14]));
% Sample time stays inherited: a root port's parameters are evaluated even
% when its choice is inactive, so fsw_period lives in the choice's hold
% Interpolation off: a replayed command holds between its rows, as the
% flight software's does
set_param(inp, 'PortDimensions', '4', 'OutDataTypeStr', 'double', 'SignalType', 'real', ...
    'Interpolate', 'off');
add_block('simulink/Sinks/Out1', [m '/x_fsw'], 'Position', blockPos(m, 'GNC', [220 -60], [30 14]));
gp = get_param(gnc, 'PortHandles');
ip = get_param(inp, 'PortHandles');
op = get_param([m '/x_fsw'], 'PortHandles');
add_line(m, ip.Outport(1), gp.Inport(end), 'autorouting', 'smart');
add_line(m, gp.Outport(end), op.Inport(1), 'autorouting', 'smart');

save_system(m);
fprintf('%s: added CONTROLLER=fsw (Controller_FSW, root u_fsw / x_fsw)\n', m);
close_system(m, 0);
end

function p = blockPos(sys, ref, offset, wh)
% position offset from block ref's top-left corner; wh = [width height]
% (default: ref's own size)
r = get_param([sys '/' ref], 'Position');
if nargin < 4, wh = r(3:4) - r(1:2); end
p = [r(1:2) + offset, r(1:2) + offset + wh];
end

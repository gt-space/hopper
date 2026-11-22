function addFswSensorPorts()
%ADDFSWSENSORPORTS  CONTROLLER=fsw: the Sensors block's output on root ports.
%   Gives the flight-software plant (fsw_bridge/build_plant_c) the shared
%   sensor model, so a SIL can fly the real navigation filter on hopper's
%   own sensors instead of truth:
%   - Navigation_Sensors/Sensors/Sensors_On also runs with CONTROLLER=fsw
%     (before: only with FDIR on or NAVIGATION not none). It holds the
%     sensor models and, with FAULT_INJECTION on, the fault gates
%     (src/sensors/sensor_frame.m), so its output is after the gates.
%   - Its meas (26 x 1), valid (9 x 1) and stamp (9 x 1), laid out as
%     src/sensors/sensor_layout.m, go through a variant subsystem
%     FSW_Sensor_Ports (active only with CONTROLLER=fsw, conditions
%     propagated) to Navigation_Sensors outports and on to root Outports
%     sens_meas, sens_valid and sens_stamp. With any other controller the
%     ports drop out of the compiled model, as x_fsw does.
%   Nothing feeds back from these ports: x_fsw and the flight are unchanged.
%   Rerunning is harmless.

m = 'hopper_6dof';
load_system(m);
ns = [m '/Navigation_Sensors'];
sp = [ns '/FSW_Sensor_Ports'];
if getSimulinkBlockHandle(sp) > 0
    fprintf('%s: FSW_Sensor_Ports already there\n', m);
    close_system(m, 0);
    return
end
k = find(variantAxes().CONTROLLER == "fsw");
assert(isscalar(k), 'addFswSensorPorts:axis', 'variantAxes has no single fsw choice');
fswCond = sprintf('VARIANT_CONTROLLER == %d', k);

% ---- Sensors on for the flight software ----
on = [ns '/Sensors/Sensors_On'];
off = [ns '/Sensors/Sensors_Off'];
cond = get_param(on, 'VariantControl');
assert(strcmp(cond, 'VARIANT_FDIR == 2 || VARIANT_NAVIGATION > 1'), 'addFswSensorPorts:cond', ...
    'Sensors_On''s condition is ''%s''; update this script for it.', cond);
cond = [cond ' || ' fswCond];
set_param(on, 'VariantControl', cond);
set_param(off, 'VariantControl', ['~(' cond ')']);

% ---- variant pass-through inside Navigation_Sensors ----
names = {'meas', 'valid', 'stamp'};
sensPos = get_param([ns '/Sensors'], 'Position');
pos = [sensPos(3) + 150, sensPos(4) + 80, sensPos(3) + 290, sensPos(4) + 200];
add_block('simulink/Ports & Subsystems/Variant Subsystem', sp, 'Position', pos);
% the library block comes with sample choices; start empty
kids = find_system(sp, 'SearchDepth', 1, 'MatchFilter', @Simulink.match.allVariants);
for i = 2:numel(kids), delete_block(kids{i}); end
ch = [sp '/FSW_Sensor_Ports_On'];
add_block('simulink/Ports & Subsystems/Subsystem', ch, 'Position', [150 40 290 160]);
delete_line(ch, 'In1/1', 'Out1/1');
delete_block([ch '/In1']);
delete_block([ch '/Out1']);
for i = 1:3
    y = 40 * i;
    add_block('simulink/Sources/In1', [sp '/' names{i}], 'Position', [40 y 70 y + 14]);
    add_block('simulink/Sinks/Out1', [sp '/sens_' names{i}], 'Position', [380 y 410 y + 14]);
    add_block('simulink/Sources/In1', [ch '/' names{i}], 'Position', [40 y 70 y + 14]);
    add_block('simulink/Sinks/Out1', [ch '/sens_' names{i}], 'Position', [200 y 230 y + 14]);
    add_line(ch, [names{i} '/1'], ['sens_' names{i} '/1']);
end
set_param(ch, 'VariantControl', fswCond);
set_param(sp, 'VariantControlMode', 'expression', 'PropagateVariantConditions', 'on', ...
    'AllowZeroVariantControls', 'on');
set_param(sp, 'Description', ['CONTROLLER=fsw only: the sensor model''s meas/valid/stamp out to ' ...
    'root sens_* ports for the flight-software plant. Built by tools/refactor/addFswSensorPorts.m.']);

% ---- Sensors -> FSW_Sensor_Ports -> Navigation_Sensors outports ----
for i = 1:3
    add_line(ns, sprintf('Sensors/%d', i), sprintf('FSW_Sensor_Ports/%d', i), 'autorouting', 'smart');
    o = [ns '/sens_' names{i}];
    add_block('simulink/Sinks/Out1', o, 'Position', [pos(3) + 120, pos(2) + 40 * i - 30, pos(3) + 150, pos(2) + 40 * i - 16]);
    add_line(ns, sprintf('FSW_Sensor_Ports/%d', i), ['sens_' names{i} '/1'], 'autorouting', 'smart');
end

% ---- root ports ----
nph = get_param(ns, 'PortHandles');
nsPos = get_param(ns, 'Position');
for i = 1:3
    o = [m '/sens_' names{i}];
    y = nsPos(4) + 30 + 40 * i;
    add_block('simulink/Sinks/Out1', o, 'Position', [nsPos(3) + 200, y, nsPos(3) + 230, y + 14]);
    oph = get_param(o, 'PortHandles');
    add_line(m, nph.Outport(end - 3 + i), oph.Inport(1), 'autorouting', 'smart');
end

save_system(m);
fprintf('%s: Sensors on for CONTROLLER=fsw; root sens_meas / sens_valid / sens_stamp\n', m);
close_system(m, 0);
end

function gateTouchdownStop()
%GATETOUCHDOWNSTOP  CONTROLLER=fsw with GROUND=strut runs on past touchdown.
%   Vehicle_Dynamics/Touchdown_Stop (subsystems/flight_logic/Touchdown_Stop.slx)
%   ends every run at the first descending ground contact. The
%   flight-software plant (fsw_bridge/build_plant_c 'Ground', 'strut')
%   needs the vehicle to come to rest on its struts instead, so the flight
%   software's own touchdown detection can be tested; the flight software
%   (or the replay's command file, or StopTime) then ends the run.
%
%   The block becomes a Variant Subsystem Touchdown_Stop with one choice,
%   Stop (the same subsystem reference), active unless
%   CONTROLLER=fsw and GROUND=strut. With no ground (GROUND=none) the stop
%   stays on for every controller: without struts the vehicle would fall
%   through the pad. Every other configuration compiles exactly as before.
%   Rerunning is harmless.

m = 'hopper_6dof';
load_system(m);
vd = [m '/Vehicle_Dynamics'];
blk = [vd '/Touchdown_Stop'];
if strcmp(get_param(blk, 'Variant'), 'on')
    fprintf('%s: Touchdown_Stop already gated\n', m);
    close_system(m, 0);
    return
end
assert(strcmp(get_param(blk, 'ReferencedSubsystem'), 'Touchdown_Stop'), 'gateTouchdownStop:ref', ...
    'Touchdown_Stop is not the subsystem reference this script expects');
A = variantAxes();
kf = find(A.CONTROLLER == "fsw");
kg = find(A.GROUND == "strut");
assert(isscalar(kf) && isscalar(kg), 'gateTouchdownStop:axis', 'variantAxes has no single fsw / strut choice');
cond = sprintf('~(VARIANT_CONTROLLER == %d && VARIANT_GROUND == %d)', kf, kg);

% ---- record the wiring, then replace the block ----
pos = get_param(blk, 'Position');
ph = get_param(blk, 'PortHandles');
src = zeros(1, numel(ph.Inport));
for i = 1:numel(ph.Inport)
    l = get_param(ph.Inport(i), 'Line');
    src(i) = get_param(l, 'SrcPortHandle');
    delete_line(l);
end
delete_block(blk);

add_block('simulink/Ports & Subsystems/Variant Subsystem', blk, 'Position', pos);
kids = find_system(blk, 'SearchDepth', 1, 'MatchFilter', @Simulink.match.allVariants);
for i = 2:numel(kids), delete_block(kids{i}); end
names = {'Z', 'W'};
for i = 1:2
    y = 40 * i;
    add_block('simulink/Sources/In1', [blk '/' names{i}], 'Position', [40 y 70 y + 14]);
end
ch = [blk '/Stop'];
add_block('simulink/Ports & Subsystems/Subsystem Reference', ch, 'Position', [150 30 290 110]);
set_param(ch, 'ReferencedSubsystem', 'Touchdown_Stop');
set_param(ch, 'VariantControl', cond);
set_param(blk, 'VariantControlMode', 'expression', 'PropagateVariantConditions', 'off', ...
    'AllowZeroVariantControls', 'on');
set_param(blk, 'Description', ['Stops the run at the first descending ground contact, except with ' ...
    'CONTROLLER=fsw and GROUND=strut (the flight software''s plant rests on its struts). ' ...
    'Built by tools/refactor/gateTouchdownStop.m.']);

nph = get_param(blk, 'PortHandles');
for i = 1:numel(src)
    add_line(vd, src(i), nph.Inport(i), 'autorouting', 'smart');
end

save_system(m);
fprintf('%s: Touchdown_Stop off for CONTROLLER=fsw with GROUND=strut\n', m);
close_system(m, 0);
end

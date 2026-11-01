function parameterizeLandingCutoff()
%PARAMETERIZELANDINGCUTOFF  Landing_Engine_Cutoff height from IN.mission.
%   subsystems/flight_logic/Landing_Engine_Cutoff.slx cut the engine below
%   a hard-coded 5 cm while descending (If: u1 > -0.05 & u2 > 0, u1 = NED
%   z). The height is now IN.mission.engine_cutoff_alt, added to z before
%   the If (u1 > 0 & u2 > 0), so a vehicle that can throttle down to the
%   ground (the CTV's EDF) can cut later. Rerunning is harmless.

f = 'Landing_Engine_Cutoff';
load_system(f);
ifb = [f '/If'];
if getSimulinkBlockHandle([f '/Cutoff_Altitude']) > 0
    fprintf('%s: already parameterized\n', f); close_system(f, 0); return
end
ph = get_param(ifb, 'PortHandles');
l = get_param(ph.Inport(1), 'Line');
src = get_param(l, 'SrcPortHandle');
delete_line(f, src, ph.Inport(1));
pos = get_param(ifb, 'Position');
add_block('built-in/Constant', [f '/Cutoff_Altitude'], 'Value', 'IN.mission.engine_cutoff_alt', ...
    'Position', [pos(1) - 160, pos(2) - 40, pos(1) - 100, pos(2) - 20]);
add_block('built-in/Sum', [f '/Height_Above_Cutoff'], 'Inputs', '++', 'IconShape', 'round', ...
    'Position', [pos(1) - 60, pos(2) + 5, pos(1) - 40, pos(2) + 25]);
s = get_param([f '/Height_Above_Cutoff'], 'PortHandles');
c = get_param([f '/Cutoff_Altitude'], 'PortHandles');
add_line(f, src, s.Inport(1), 'autorouting', 'smart');
add_line(f, c.Outport, s.Inport(2), 'autorouting', 'smart');
add_line(f, s.Outport, ph.Inport(1), 'autorouting', 'smart');
set_param(ifb, 'IfExpression', 'u1 > 0 & u2 > 0');   % z + cutoff > 0: below the cutoff height (NED)
save_system(f);
close_system(f);
fprintf('%s: cutoff height now IN.mission.engine_cutoff_alt\n', f);
end

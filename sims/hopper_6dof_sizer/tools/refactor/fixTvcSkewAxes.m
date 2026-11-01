function fixTvcSkewAxes()
%FIXTVCSKEWAXES  Engine-CG swing uses the same gimbal axes as the thrust.
%   subsystems/mass_properties/Mass_Properties.slx, TVC_Skew_Matrix: the
%   small-rotation matrix [theta]x that swings the engine CG about the
%   gimbal pivot (theta x r, r = engine CG - pivot). It took its inputs as
%   theta = [0; delta_y; delta_z] with delta_y wired to the TVC *yaw*
%   angle and delta_z to pitch. The thrust model (Actuators_TVC_RCS) tilts
%   the thrust by pitch about body y and yaw about z. So the engine CG
%   swung off its own thrust line, and cross(arm, F) made roll torque a
%   pivoted engine can't (about 0.04 N m on the CTV) and mixed the pitch
%   and yaw moments. Present since v2.
%   The block now names its inputs pitch and yaw, with theta = [0; pitch;
%   yaw], and this script checks against the wiring that the input named
%   pitch is fed by the TVC Pitch Angle inport (Simulink keeps a line on
%   its named input when the inputs are reordered, so port order alone
%   can't be trusted). Rerunning is harmless.

f = 'Mass_Properties';
load_system(f);
blk = [f '/TVC_Skew_Matrix'];
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
if ~contains(ch.Script, 'function rot_mat = fcn(pitch, yaw)')
    ch.Script = strjoin({
        'function rot_mat = fcn(pitch, yaw)'
        '% [theta]x for the gimbal rotation theta = [0; pitch; yaw]: pitch about'
        '% body y and yaw about body z, as Actuators_TVC_RCS tilts the thrust,'
        '% so the engine CG swings along its own thrust line'
        '% (tools/refactor/fixTvcSkewAxes.m)'
        ''
        'rot_mat = [0 -yaw pitch;'
        '           yaw 0 0;'
        '           -pitch 0 0];'
        'end'}, newline);
    changed = true;
else
    changed = false;
end
% each named input must be fed by the inport of the same angle
ph = get_param(blk, 'PortHandles');
for d = ch.find('-isa', 'Stateflow.Data', 'Scope', 'Input')'
    src = get_param(get_param(get_param(ph.Inport(d.Port), 'Line'), 'SrcPortHandle'), 'Parent');
    want = struct('pitch', 'TVC Pitch Angle', 'yaw', 'TVC Yaw Angle');
    got = get_param(src, 'Name');
    if ~strcmp(got, want.(d.Name))
        error('fixTvcSkewAxes:wiring', '%s input %s is fed by %s, not %s.', blk, d.Name, got, want.(d.Name));
    end
end
if changed
    save_system(f);
    fprintf('%s: theta = [0; pitch; yaw], inputs checked against the wiring\n', blk);
else
    fprintf('%s: already fixed (wiring checked)\n', blk);
end
close_system(f, 0);
end

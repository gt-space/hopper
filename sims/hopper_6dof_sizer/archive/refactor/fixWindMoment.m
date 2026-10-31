function fixWindMoment(model)
%FIXWINDMOMENT  Wind_Forces: the wind moment as r_cp x F.
%   fixWindMoment()        fix subsystems/environment/Wind_Forces.slx and
%                          wire hopper_6dof_NED_v2
%   fixWindMoment(model)   wire another model that references the file
%
%   Wind_Forces put out the wind moment as F .* (cp - cg), element by
%   element, in inertial axes: a north force gave a moment about north (it
%   belongs about east), and the east force used cp itself as its arm. The
%   moment is now cross(r, F) with r the centre of pressure from the CG,
%   cp - cg_x along the body x axis, rotated to inertial axes with the
%   transpose of the DCM. The model still rotates force and moment into
%   body axes after the block (MatrixMultiply1/2), so the force path is
%   unchanged. The block gains a third input, 'Inertial to Body' (the
%   6DOF block's DCM). Rerunning is harmless.

if nargin < 1, model = 'hopper_6dof_NED_v2'; end
f = 'Wind_Forces';
load_system(f);
if getSimulinkBlockHandle([f '/Wind_Moment']) < 0
    % replace the element-wise moment with r x F
    for b = {'Product3', 'Subtract2'}
        lh = get_param([f '/' b{1}], 'LineHandles');
        ph = get_param([f '/' b{1}], 'PortHandles');
        for k = 1:numel(lh.Inport)
            if lh.Inport(k) > 0, delete_line(f, get_param(lh.Inport(k), 'SrcPortHandle'), ph.Inport(k)); end
        end
        for h = lh.Outport(lh.Outport > 0), delete_line(h); end
        delete_block([f '/' b{1}]);
    end
    pos = get_param([f '/Mxyz_wind'], 'Position');
    add_block('built-in/Inport', [f '/Inertial to Body'], 'Port', '3', 'Position', pos + [-420 60 -420 60]);
    fn = [f '/Wind_Moment'];
    add_block('simulink/User-Defined Functions/MATLAB Function', fn, 'Position', pos + [-220 -30 -120 70]);
    ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', fn);
    ch.Script = strjoin({
        'function M = fcn(F, DCM, cg, cp)'
        '% Wind moment about the CG in inertial axes (the model rotates it to'
        '% body axes with the DCM, like the force): r x F, with the centre of'
        '% pressure cp metres from the base along the body x axis'
        'r = DCM'' * [cp - cg(1); 0; 0];'
        'M = cross(r, F(:));'
        'end'}, newline);
    add_line(f, 'Mux/1', 'Wind_Moment/1', 'autorouting', 'smart');
    add_line(f, 'Inertial to Body/1', 'Wind_Moment/2', 'autorouting', 'smart');
    add_line(f, 'cg/1', 'Wind_Moment/3', 'autorouting', 'smart');
    add_line(f, 'cp/1', 'Wind_Moment/4', 'autorouting', 'smart');
    add_line(f, 'Wind_Moment/1', 'Mxyz_wind/1', 'autorouting', 'smart');
    save_system(f);
    fprintf('%s: moment is now r x F\n', f);
end

% the model feeds the DCM it already uses to rotate the wind into body axes
load_system(model);
wf = get_param([model '/Wind_Forces'], 'PortHandles');
if get_param(wf.Inport(3), 'Line') < 0
    six = get_param([model '/6DOF (Quaternion)'], 'PortHandles');
    add_line(model, six.Outport(5), wf.Inport(3), 'autorouting', 'smart');
    save_system(model);
    fprintf('%s: DCM wired into Wind_Forces\n', model);
end
end

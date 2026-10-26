function converted = convertToSubsystemReferences(model)
%CONVERTTOSUBSYSTEMREFERENCES  Move v2's shared blocks into subsystem files.
%   converted = convertToSubsystemReferences()       hopper_6dof_NED_v2
%   converted = convertToSubsystemReferences(model)
%
%   Each block listed below becomes a Subsystem Reference to a .slx file
%   under subsystems/. Every model that references the file uses the same
%   copy, so a fix made in the file reaches all of them. Blocks that are
%   already references are skipped, so the script can be rerun as the list
%   grows (stage 5, step 2 onward). Contents and wiring don't change.
%   Returns {block, file} for every block converted.
%
%   Several blocks can share one file: the six slosh pendulums (ox/fuel x
%   axial/pitch/yaw) all reference slosh/Slosh_Axis.slx. Their copies
%   differed only in port and variable names, so the first one is given
%   generic names (accel, damping, pendulum_l, ...) before it becomes the
%   file. Nested blocks come before the block that contains them.

if nargin < 1, model = 'hopper_6dof_NED_v2'; end

% {block (relative to the model), file under subsystems/}
map = {
    'Propellant_Mass_Integration',      'propulsion/Propellant_Mass_Integration.slx'
    'Thrust_to_Mdot_LUT',               'propulsion/Thrust_to_Mdot_LUT.slx'
    'Propellant_Depletion_Cutoff',      'flight_logic/Propellant_Depletion_Cutoff.slx'
    'Landing_Engine_Cutoff',            'flight_logic/Landing_Engine_Cutoff.slx'
    'Touchdown_Stop',                   'flight_logic/Touchdown_Stop.slx'
    'Slosh/Oxidizer Axial Force',       'slosh/Slosh_Axis.slx'
    'Slosh/Oxidizer Pitch Force',       'slosh/Slosh_Axis.slx'
    'Slosh/Oxidizer Yaw Force',         'slosh/Slosh_Axis.slx'
    'Slosh/Fuel Axial Force',           'slosh/Slosh_Axis.slx'
    'Slosh/Fuel Pitch Force',           'slosh/Slosh_Axis.slx'
    'Slosh/Fuel Yaw Force',             'slosh/Slosh_Axis.slx'
    'Slosh',                            'slosh/Slosh.slx'
    'Wind_Forces',                      'environment/Wind_Forces.slx'
    'Actuators_TVC_RCS/EDF_Dynamics',   'actuators/EDF_Dynamics.slx'
    'Actuators_TVC_RCS',                'actuators/Actuators_TVC_RCS.slx'
    'Mass_Properties',                  'mass_properties/Mass_Properties.slx'
    };

load_system(model);
proj = currentProject();
converted = {};
for k = 1:size(map, 1)
    blk = [model '/' map{k, 1}];
    if getSimulinkBlockHandle(blk) < 0
        fprintf('%s: no %s, skipped\n', model, map{k, 1}); continue
    end
    if ~isempty(get_param(blk, 'ReferencedSubsystem'))
        fprintf('%s: %s already references %s\n', model, map{k, 1}, get_param(blk, 'ReferencedSubsystem')); continue
    end
    file = hopperRoot('subsystems', map{k, 2});
    [folder, name] = fileparts(file);
    if ~isfolder(folder)
        mkdir(folder);
        addPath(proj, folder);
    end
    if isfile(file)
        % Another block or model already created the shared file: point at it
        pointAtFile(blk, name);
    else
        if strcmp(name, 'Slosh_Axis')
            genericSloshAxisNames(blk);
        end
        Simulink.SubsystemReference.convertSubsystemToSubsystemReference(blk, file);
        addFile(proj, file);
    end
    converted(end+1, :) = {map{k, 1}, map{k, 2}}; %#ok<AGROW>
    fprintf('%s: %s -> subsystems/%s\n', model, map{k, 1}, map{k, 2});
end
if ~isempty(converted)
    save_system(model);
end
end

function pointAtFile(blk, name)
% Replace blk's contents with a reference to subsystem file name, keeping
% its wiring. Emptying the block removes its ports, so the lines are
% recorded first and redrawn to the same port numbers afterwards.
sys = get_param(blk, 'Parent');
ph = get_param(blk, 'PortHandles');
ins = {}; outs = {};
for p = 1:numel(ph.Inport)
    l = get_param(ph.Inport(p), 'Line');
    if l > 0
        ins(end+1, :) = {get_param(l, 'SrcPortHandle'), p}; %#ok<AGROW>
        delete_line(sys, ins{end, 1}, ph.Inport(p));
    end
end
for p = 1:numel(ph.Outport)
    l = get_param(ph.Outport(p), 'Line');
    if l > 0
        outs(end+1, :) = {p, get_param(l, 'DstPortHandle'), get_param(l, 'Name')}; %#ok<AGROW>
        delete_line(l);
    end
end
Simulink.SubSystem.deleteContents(blk);
set_param(blk, 'ReferencedSubsystem', name);
ph = get_param(blk, 'PortHandles');
for r = 1:size(ins, 1)
    add_line(sys, ins{r, 1}, ph.Inport(ins{r, 2}), 'autorouting', 'smart');
end
for r = 1:size(outs, 1)
    dst = outs{r, 2};
    for d = dst(dst > 0)'
        l = add_line(sys, ph.Outport(outs{r, 1}), d, 'autorouting', 'smart');
        set_param(l, 'Name', outs{r, 3});
    end
end
end

function genericSloshAxisNames(blk)
% ox_axial_accel -> accel, z_omega_dot -> omega_dot, fu_slosh_mass_y ->
% bob_mass, ... in the axis block's ports, MATLAB Function code and line
% labels. Only names change; the MATLAB Function ports follow their code.
generic = @(s) regexprep(s, { ...
    '\<(ox|fu)_slosh_mass_[xyz]\>', '\<(ox|fu)_slosh_torque_[xyz]\>', '\<(ox|fu)_slosh_[xyz]\>', ...
    '\<(ox|fu)_(axial|pitch|yaw)_f\>', '\<(ox|fu)_(axial|lateral|pitch|yaw)_', '\<[xyz]_(omega_dot|slosh_torque)\>'}, ...
    {'bob_mass', 'hinge_torque', 'bob_position', 'force', '', '$1'});
for b = [find_system(blk, 'SearchDepth', 1, 'BlockType', 'Inport'); find_system(blk, 'SearchDepth', 1, 'BlockType', 'Outport')]'
    set_param(b{1}, 'Name', generic(get_param(b{1}, 'Name')));
end
charts = sfroot().find('-isa', 'Stateflow.EMChart');
for c = charts(startsWith({charts.Path}, [blk '/']))'
    c.Script = generic(c.Script);
end
for l = find_system(blk, 'FindAll', 'on', 'Type', 'line')'
    n = get_param(l, 'Name');
    if ~isempty(n), set_param(l, 'Name', generic(n)); end
end
end

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
%   file. Nested blocks come before the block that contains them. The fuel
%   and ox throttle valve controllers share propulsion/valves/
%   Throttle_Valve_Controller.slx, with the one setting that differs
%   (ramp_slope) as a mask parameter set per instance. MATLAB Function
%   blocks (LQR_Gain_Scheduler) are first wrapped in a plain subsystem.
%   Entries with a model name apply only to that model.

if nargin < 1, model = 'hopper_6dof_NED_v2'; end

% {block (relative to the model), file under subsystems/, only in this model}
% '' = v2 and the models built from it. Other feature copies (such as
% hopper_6dof_NED_v2_ground) hold older versions of some blocks, so only
% the entries that name them apply there.
fluid = 'hopper_6dof_NED_v2_fluid_integration';
fromV2 = {'hopper_6dof_NED_v2', fluid};
map = {
    'Propellant_Mass_Integration',      'propulsion/Propellant_Mass_Integration.slx',  ''
    'Thrust_to_Mdot_LUT',               'propulsion/Thrust_to_Mdot_LUT.slx',           ''
    'Propellant_Depletion_Cutoff',      'flight_logic/Propellant_Depletion_Cutoff.slx', ''
    'Landing_Engine_Cutoff',            'flight_logic/Landing_Engine_Cutoff.slx',      ''
    'Touchdown_Stop',                   'flight_logic/Touchdown_Stop.slx',             ''
    'Slosh/Oxidizer Axial Force',       'slosh/Slosh_Axis.slx',                        ''
    'Slosh/Oxidizer Pitch Force',       'slosh/Slosh_Axis.slx',                        ''
    'Slosh/Oxidizer Yaw Force',         'slosh/Slosh_Axis.slx',                        ''
    'Slosh/Fuel Axial Force',           'slosh/Slosh_Axis.slx',                        ''
    'Slosh/Fuel Pitch Force',           'slosh/Slosh_Axis.slx',                        ''
    'Slosh/Fuel Yaw Force',             'slosh/Slosh_Axis.slx',                        ''
    'Slosh',                            'slosh/Slosh.slx',                             ''
    'Wind_Forces',                      'environment/Wind_Forces.slx',                 ''
    'Actuators_TVC_RCS/EDF_Dynamics',   'actuators/EDF_Dynamics.slx',                  ''
    'Actuators_TVC_RCS',                'actuators/Actuators_TVC_RCS.slx',             ''
    'Mass_Properties',                  'mass_properties/Mass_Properties.slx',         ''
    'LQR_Gain_Scheduler',               'gnc/LQR_Gain_Scheduler.slx',                  ''
    % v2's commented-out leg model (scopes only): four identical legs
    'Ground_Contact_Unwired/Leg_1_Force', 'ground/Leg_Force_Unwired.slx',              ''
    'Ground_Contact_Unwired/Leg_2_Force', 'ground/Leg_Force_Unwired.slx',              ''
    'Ground_Contact_Unwired/Leg_3_Force', 'ground/Leg_Force_Unwired.slx',              ''
    'Ground_Contact_Unwired/Leg_4_Force', 'ground/Leg_Force_Unwired.slx',              ''
    % #27's strut ground contact: 'Ground' in hopper_6dof_NED_v2_ground,
    % copied into the fluid model as Ground_Contact
    'Ground',                           'ground/Ground_Contact.slx',                   'hopper_6dof_NED_v2_ground'
    'Ground_Contact',                   'ground/Ground_Contact.slx',                   fluid
    % #26's fluid network exists live only in the fluid model; v2 keeps an
    % older, commented-out version with different valve controllers
    'Fluid_Network/Fuel Throttle Valve Controller', 'propulsion/valves/Throttle_Valve_Controller.slx', fluid
    'Fluid_Network/Ox Throttle Valve Controller',   'propulsion/valves/Throttle_Valve_Controller.slx', fluid
    };
% Mask parameter values for instances of a shared file that differ per use
instanceParams = containers.Map();
instanceParams('Fluid_Network/Fuel Throttle Valve Controller') = {'ramp_slope', '60'};
instanceParams('Fluid_Network/Ox Throttle Valve Controller')   = {'ramp_slope', '48'};

load_system(model);
proj = currentProject();
converted = {};
for k = 1:size(map, 1)
    if isempty(map{k, 3})
        if ~ismember(model, fromV2), continue; end
    elseif ~strcmp(map{k, 3}, model)
        continue
    end
    blk = [model '/' map{k, 1}];
    if getSimulinkBlockHandle(blk) < 0
        fprintf('%s: no %s, skipped\n', model, map{k, 1}); continue
    end
    if strcmp(get_param(blk, 'SFBlockType'), 'MATLAB Function')
        blk = wrapInSubsystem(blk);
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
        switch name
            case 'Slosh_Axis', genericSloshAxisNames(blk);
            case 'Throttle_Valve_Controller', commonValveController(blk);
        end
        Simulink.SubsystemReference.convertSubsystemToSubsystemReference(blk, file);
        addFile(proj, file);
        if strcmp(name, 'Throttle_Valve_Controller')
            maskValveController(name);
        end
    end
    if isKey(instanceParams, map{k, 1})
        nv = instanceParams(map{k, 1});
        set_param(blk, nv{:});
    end
    converted(end+1, :) = {map{k, 1}, map{k, 2}}; %#ok<AGROW>
    fprintf('%s: %s -> subsystems/%s\n', model, map{k, 1}, map{k, 2});
end
if ~isempty(converted)
    save_system(model);
end
end

function blk = wrapInSubsystem(blk)
% A MATLAB Function block can't be a Subsystem Reference itself. Put it in
% a plain subsystem that takes over its name and position; the function
% inside keeps its own name (LQR_Gain_Scheduler holds Gain_Selector).
sys = get_param(blk, 'Parent');
name = get_param(blk, 'Name');
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
inner = regexp(ch.Script, 'function\s+[^=]*=\s*(\w+)', 'tokens', 'once');
d = ch.find('-isa', 'Stateflow.Data');
argIn = portOrder(d(strcmp({d.Scope}, 'Input')));
argOut = portOrder(d(strcmp({d.Scope}, 'Output')));
pos = get_param(blk, 'Position');
before = find_system(sys, 'SearchDepth', 1, 'BlockType', 'SubSystem');
Simulink.BlockDiagram.createSubsystem(get_param(blk, 'Handle'));
wrapper = setdiff(find_system(sys, 'SearchDepth', 1, 'BlockType', 'SubSystem'), before);
wrapper = wrapper{1};
set_param([wrapper '/' name], 'Name', inner{1});   % the function, now inside
set_param(wrapper, 'Name', name, 'Position', pos);
blk = [sys '/' name];
% wrapper ports take the function's argument names
ins = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Inport');
outs = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Outport');
for p = 1:numel(ins)
    set_param(ins{p}, 'Name', argIn{str2double(get_param(ins{p}, 'Port'))});
end
for p = 1:numel(outs)
    set_param(outs{p}, 'Name', argOut{str2double(get_param(outs{p}, 'Port'))});
end
end

function commonValveController(blk)
% #26's fuel and ox throttle valve controllers are the same design. They
% differ in the start-up ramp slope (60 vs 48, now the ramp_slope mask
% parameter), the ramp hand-off switch (fuel u2 > 0.8 s, ox u2 >= 0.8 s;
% the fuel version is kept), one setpoint port name (P_set vs pset) and
% debug logging. The To Workspace blocks (theta_log, Pout_log,
% feed_forward) are removed: nothing reads them from this model, and two
% instances would write the same variables.
ramp = find_system(blk, 'LookUnderMasks', 'all', 'Name', 'Ramp');
set_param(ramp{1}, 'slope', 'ramp_slope');
for tw = find_system(blk, 'LookUnderMasks', 'all', 'BlockType', 'ToWorkspace')'
    lh = get_param(tw{1}, 'LineHandles');
    ph = get_param(tw{1}, 'PortHandles');
    if lh.Inport > 0   % remove only this branch of the logged signal
        delete_line(get_param(tw{1}, 'Parent'), get_param(lh.Inport, 'SrcPortHandle'), ph.Inport);
    end
    delete_block(tw{1});
end
end

function maskValveController(file)
load_system(file);
mk = Simulink.Mask.create(file);
mk.addParameter('Name', 'ramp_slope', 'Prompt', 'Start-up ramp slope', 'Type', 'edit', 'Value', '60', 'Evaluate', 'on', 'Tunable', 'on');
mk.Description = 'Throttle valve angle controller (PR #26). Shared by the fuel and ox valves; ramp_slope is set per instance.';
save_system(file);
end

function names = portOrder(data)
names = cell(1, numel(data));
for k = 1:numel(data)
    names{data(k).Port} = data(k).Name;
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

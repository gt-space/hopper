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

if nargin < 1, model = 'hopper_6dof_NED_v2'; end

% {block (relative to the model), file under subsystems/}
map = {
    'Propellant_Mass_Integration',  'propulsion/Propellant_Mass_Integration.slx'
    'Thrust_to_Mdot_LUT',           'propulsion/Thrust_to_Mdot_LUT.slx'
    'Propellant_Depletion_Cutoff',  'flight_logic/Propellant_Depletion_Cutoff.slx'
    'Landing_Engine_Cutoff',        'flight_logic/Landing_Engine_Cutoff.slx'
    'Touchdown_Stop',               'flight_logic/Touchdown_Stop.slx'
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
    folder = fileparts(file);
    if ~isfolder(folder)
        mkdir(folder);
        addPath(proj, folder);
    end
    if isfile(file)
        % Another model already created the shared file: point at it
        Simulink.SubSystem.deleteContents(blk);
        [~, name] = fileparts(file);
        set_param(blk, 'ReferencedSubsystem', name);
    else
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

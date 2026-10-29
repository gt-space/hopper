function parameterizeVehicleGeometry()
%PARAMETERIZEVEHICLEGEOMETRY  Shared subsystems read geometry from the vehicle config.
%   Tank positions were typed into the Slosh and Mass_Properties blocks as
%   [-(1.0668 - 18.91*0.0254/2) 0 0] (the hopper's tank bottoms), the fan
%   moment arm was the hopper's tank radius, and the ground contact read
%   the sizer's leg stiffness. They now read IN.vehicle (the config,
%   vehicles/<vehicle>/<vehicle>_config.m, via vehicleInputs), so another
%   vehicle gets its own geometry. The expressions evaluate to the same
%   numbers for the hopper. Rerunning is harmless.

edits = {
    % file                 old value                              new value
    'Slosh',             '[-(1.0668 - 18.91*0.0254/2) 0 0]',    '[-IN.vehicle.tanks.fu.floor_x 0 0]'
    'Slosh',             '[-(1.854 - 29.18*0.0254/2) 0 0]',     '[-IN.vehicle.tanks.ox.floor_x 0 0]'
    'Slosh',             'TANKS.singular.radius',               'IN.vehicle.tanks.radius'
    'Mass_Properties',   '[(1.0668 - 18.91*0.0254/2) 0 0]',     '[IN.vehicle.tanks.fu.floor_x 0 0]'
    'Mass_Properties',   '[(1.854 - 29.18*0.0254/2) 0 0]',      '[IN.vehicle.tanks.ox.floor_x 0 0]'
    'EDF_Dynamics',      'TANKS.singular.radius',               'IN.vehicle.rcs.arm'
    'FanModule_2025a',   'TANKS.singular.radius',               'IN.vehicle.rcs.arm'
    'Ground_Contact',    'STRUCT.legs.k',                       'IN.vehicle.legs.stiffness'
    'Ground_Contact',    'r_leg',                               'IN.vehicle.legs.r_foot'
    };
for f = unique(edits(:, 1), 'stable')'
    load_system(f{1});
    n = 0;
    for k = find(strcmp(edits(:, 1), f{1}))'
        b = find_system(f{1}, 'LookUnderMasks', 'all', 'FollowLinks', 'off', 'BlockType', 'Constant', 'Value', edits{k, 2});
        for i = 1:numel(b)
            set_param(b{i}, 'Value', edits{k, 3});
            n = n + 1;
        end
        % mask or block parameters other than a Constant's Value
        if isempty(b)
            n = n + replaceInDialogs(f{1}, edits{k, 2}, edits{k, 3});
        end
    end
    if n > 0
        save_system(f{1});
    end
    close_system(f{1}, 0);
    fprintf('%s: %d parameter(s) now read IN.vehicle\n', f{1}, n);
end
end

function n = replaceInDialogs(sys, old, new)
% Any block dialog parameter whose whole value is old
n = 0;
blocks = find_system(sys, 'LookUnderMasks', 'all', 'FollowLinks', 'off', 'Type', 'Block');
for i = 1:numel(blocks)
    dp = get_param(blocks{i}, 'DialogParameters');
    if ~isstruct(dp), continue; end
    for p = fieldnames(dp)'
        try
            v = get_param(blocks{i}, p{1});
        catch
            continue
        end
        if ischar(v) && strcmp(strtrim(v), old)
            set_param(blocks{i}, p{1}, new);
            n = n + 1;
        end
    end
end
end

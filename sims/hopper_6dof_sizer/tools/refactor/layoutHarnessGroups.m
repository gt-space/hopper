function layoutHarnessGroups(model)
%LAYOUTHARNESSGROUPS  Lay out the regrouped harness: groups left to right, each arranged inside.
%   layoutHarnessGroups()        models/hopper_6dof.slx (saves it)
%   layoutHarnessGroups(model)   another loaded copy (saves it)
%
%   Top level, along the signal flow: Environment, Navigation_Sensors,
%   GNC, Actuation, Vehicle_Dynamics in a row, with Mass_Propulsion and
%   Logging on a row below (the plant's feedback and the outputs). Inside
%   each group Simulink.BlockDiagram.arrangeSystem lays the blocks out; a
%   group it fails on keeps its layout and is reported. Moving blocks and
%   rerouting lines doesn't change the simulation. Rerunning is harmless.

if nargin < 1, model = 'hopper_6dof'; end
load_system(model);
row1 = {'Environment', 'Navigation_Sensors', 'GNC', 'Actuation', 'Vehicle_Dynamics'};
row2 = {'', '', 'Mass_Propulsion', '', 'Logging'};
for g = [row1, row2(~cellfun(@isempty, row2))]
    try
        Simulink.BlockDiagram.arrangeSystem([model '/' g{1}]);
    catch err
        fprintf('  %s: not auto-arranged (%s)\n', g{1}, extractBefore([err.message newline], newline));
    end
end
w = 300; gap = 220; x0 = 80; y0 = 80;
for k = 1:numel(row1)
    place([model '/' row1{k}], x0 + (k - 1) * (w + gap), y0, w);
end
y1 = y0 + max(cellfun(@(g) height([model '/' g]), row1)) + 260;
for k = find(~cellfun(@isempty, row2))
    place([model '/' row2{k}], x0 + (k - 1) * (w + gap), y1, w);
end
if isempty(find_system(model, 'SearchDepth', 2, 'BlockType', 'Outport', 'IsBusElementPort', 'on'))
    Simulink.BlockDiagram.routeLine(find_system(model, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line'));
else
    Simulink.BlockDiagram.arrangeSystem(model);   % on buses: few enough lines for a clean auto layout
end
save_system(model);
fprintf('%s: groups laid out\n', model);
end

function place(blk, x, y, w)
set_param(blk, 'Position', [x, y, x + w, y + height(blk)]);
end

function h = height(blk)
p = get_param(blk, 'Ports');
h = max(160, 28 * max(p(1:2)));
end

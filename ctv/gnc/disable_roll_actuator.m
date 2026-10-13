function disable_roll_actuator(varargin)
% DISABLE_ROLL_ACTUATOR  Remove the EDF/RCS roll actuator from a ctv model.
%
%   disable_roll_actuator
%   disable_roll_actuator('Model', 'hopper_6dof_NED_v3_ctv')
%   disable_roll_actuator('Save', false)      % dry run
%   disable_roll_actuator('Undo', true)       % put the actuator back
%
% WHY THIS IS NEEDED
% The FDIR doc says the vehicle has no roll actuator, which we expressed as
% IN.rcs.max_thrust = 0. That alone is NOT safe: inside the EDF Dynamics
% subsystem the signal path is
%
%     Cmd Torque --> Divide --> Saturation --> Transfer Fcn --> Product1 -->
%                      ^                                          ^
%     max Thrust x Radius --------------------------------------- /
%
% i.e. the commanded torque is normalised by (max_thrust * radius) and then
% de-normalised again. With max_thrust = 0 that divisor is 0, and because
% the LQR's 4th row is also zeroed the numerator is 0 too, so Divide
% computes 0/0 = NaN. The Saturation block in that path is set to
% "Commented = through" so it does not clamp the NaN, and the NaN
% propagates into Real Torque and from there into the 6DOF moments.
%
% So the actuator has to be removed structurally, not just zeroed. This
% script sets Commented = 'on' on the EDF Dynamics subsystem, which makes
% Simulink treat it as absent and ground its output at exactly zero.
%
% NOTE: 'on' is required, NOT 'through'. EDF Dynamics is one-in/one-out, so
% 'through' would pass Cmd Torque directly to Real Torque and hand the
% controller UNLIMITED roll authority -- the opposite of the intent.
%
% The .slx is never hand-edited; this goes through set_param only.
%
% The "Roll Torque" band-limited white-noise disturbance block is left
% alone, so roll is free to drift once that block is itself re-enabled (it
% ships with Commented = 'on' as well).

p = inputParser;
p.addParameter('Model', 'hopper_6dof_NED_v2_ctv', @(s) ischar(s) || isstring(s));
p.addParameter('Save', true, @islogical);
p.addParameter('Undo', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

mdl      = char(opt.Model);
ctv_dir  = fileparts(fileparts(mfilename('fullpath')));
mdl_file = fullfile(ctv_dir, [mdl '.slx']);
assert(isfile(mdl_file), 'Model not found: %s', mdl_file);

fprintf('Loading %s ...\n', mdl_file);
load_system(mdl_file);

blk = locate_rcs(mdl);
fprintf('  roll actuator subsystem: %s\n', blk);

if opt.Undo
    set_param(blk, 'Commented', 'off');
    fprintf('  Commented -> off (roll actuator RESTORED)\n');
else
    set_param(blk, 'Commented', 'on');
    fprintf('  Commented -> on  (roll actuator REMOVED, output grounded)\n');
end

if opt.Save
    save_system(mdl);
    fprintf('Saved %s\n', mdl_file);
else
    fprintf('Save skipped (Save=false); model left dirty in memory.\n');
end

end

% ========================================================================
function blk = locate_rcs(mdl)
% Find the subsystem by what it contains -- the Constant that reads
% IN.rcs.max_thrust -- rather than by its name, so a rename does not break
% this and we cannot comment out the wrong block.

cs  = find_system(mdl, 'BlockType', 'Constant');
hit = {};
for i = 1:numel(cs)
    if strcmp(strtrim(get_param(cs{i}, 'Value')), 'IN.rcs.max_thrust')
        hit{end+1} = get_param(cs{i}, 'Parent');  %#ok<AGROW>
    end
end
hit = unique(hit);

assert(~isempty(hit), ...
    ['Could not find a subsystem containing a Constant valued ' ...
     '''IN.rcs.max_thrust''. Model structure has changed; disable the ' ...
     'roll actuator by hand and re-check this script.']);
assert(isscalar(hit), ...
    'Found %d candidate roll-actuator subsystems: %s', ...
    numel(hit), strjoin(hit, ', '));

blk = hit{1};

% sanity: it must be a subsystem, not the model root, and it should hold
% the normalising divide that causes the NaN
assert(strcmp(get_param(blk, 'Type'), 'block') && ...
       strcmp(get_param(blk, 'BlockType'), 'SubSystem'), ...
    '%s is not a subsystem; refusing to comment it out.', blk);

divs = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Product', 'Inputs', '*/');
if isempty(divs)
    warning('disable_roll_actuator:noDivide', ...
        ['%s has no normalising divide any more -- the NaN this script ' ...
         'works around may already be fixed. Commenting it out anyway, ' ...
         'since the vehicle still has no roll actuator.'], blk);
end

end

function add_onboard_gfold(varargin)
% ADD_ONBOARD_GFOLD  Install the onboard G-FOLD re-planner into
%                    ctv/hopper_6dof_NED_v3_ctv.slx using the Simulink API.
%
%   add_onboard_gfold
%   add_onboard_gfold('Save', false)        % dry run, leaves the model dirty
%   add_onboard_gfold('SampleTime', '1')    % re-plan block sample time
%
% The .slx is NEVER hand-edited. Everything here goes through add_block,
% add_line, delete_line and set_param, so the change is reviewable as code
% and repeatable. Re-running the script is safe: it removes a previously
% installed block and its wiring first.
%
% WHAT IT DOES
%   Before:  Constant8 (the reference table)   --> Gain_Selector/xref
%            Constant7 (unom)                   --> Gain_Selector/unom
%
%   After:   Constant8 ----\
%            Constant7 -----> [GFOLD Replan] --> Gain_Selector/xref
%            Constant4 ----/        |         \-> Gain_Selector/unom
%            Clock, state, mass ----/         \-> gfold_flag (F-23)
%
%   The Gain_Selector MATLAB Function itself is untouched: it still receives
%   an N x 13 reference table and an N x 4 nominal control table and still
%   interpolates them by time. The replan block simply overwrites the
%   descent rows of those tables with a freshly solved G-FOLD trajectory.
%   K1flat / K2grid still come straight from the setup script.
%
% The block's body is ctv/gnc/gfold_replan.m, read in with fileread, so the
% script that lives in the model and the file in the repo cannot diverge.
%
% Defaults to hopper_6dof_NED_v2_ctv, the active CTV model. Works on v3 too
% via 'Model'; the reference-table Constant is called Target_Trajectory in
% v2 and Target_Trajectory_full in v3, and both are accepted.
%
% Run disable_roll_actuator first -- unrelated to this change, but the model
% will produce NaN without it.

p = inputParser;
p.addParameter('Model', 'hopper_6dof_NED_v2_ctv', @(s) ischar(s) || isstring(s));
p.addParameter('Save', true, @islogical);
p.addParameter('SampleTime', '1', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

mdl      = char(opt.Model);
gnc_dir  = fileparts(mfilename('fullpath'));
ctv_dir  = fileparts(gnc_dir);
mdl_file = fullfile(ctv_dir, [mdl '.slx']);
body_fil = fullfile(gnc_dir, 'gfold_replan.m');

assert(isfile(mdl_file), 'Model not found: %s', mdl_file);
assert(isfile(body_fil), 'Block body not found: %s', body_fil);

addpath(gnc_dir);

BLK   = [mdl '/GFOLD Replan'];
FLAG  = [mdl '/gfold_flag'];

fprintf('Loading %s ...\n', mdl_file);
load_system(mdl_file);

%% ---- locate the blocks we need -----------------------------------------
sel = locate(mdl);

fprintf('  Gain_Selector subsystem : %s\n', sel.gain);
fprintf('  tgrid constant          : %s\n', sel.tgrid);
fprintf('  xref constant           : %s\n', sel.xref);
fprintf('  unom constant           : %s\n', sel.unom);
fprintf('  clock                   : %s\n', sel.clock);
fprintf('  state mux               : %s\n', sel.state);
fprintf('  mass signal             : %s\n', sel.mass);

%% ---- make the script idempotent ----------------------------------------
if ~isempty(find_system(mdl, 'SearchDepth', 1, 'Name', 'GFOLD Replan'))
    fprintf('Existing GFOLD Replan block found -- removing and reinstalling.\n');
    delete_block(BLK);
    if ~isempty(find_system(mdl, 'SearchDepth', 1, 'Name', 'gfold_flag'))
        delete_block(FLAG);
    end
end

% whatever currently drives xref / unom on the Gain_Selector has to go
drop_line(sel.gain, 4);
drop_line(sel.gain, 7);

%% ---- add the re-plan block ---------------------------------------------
fprintf('Adding GFOLD Replan block ...\n');
gpos = get_param(sel.gain, 'Position');
add_block('simulink/User-Defined Functions/MATLAB Function', BLK, ...
    'Position', [gpos(1)-260, gpos(2), gpos(1)-150, gpos(4)]);

% install the body from the repo file
chart = emchart(BLK);
chart.Script = fileread(body_fil);

% the block's ports come from the script signature, so they only exist once
% the script is installed
ph = get_param(BLK, 'PortHandles');
assert(numel(ph.Inport) == 6 && numel(ph.Outport) == 3, ...
    ['GFOLD Replan has %d inputs / %d outputs, expected 6 / 3. The ' ...
     'signature in gfold_replan.m and this script disagree.'], ...
    numel(ph.Inport), numel(ph.Outport));

% low re-plan rate. gfold_onboard_step also gates itself to 1 Hz of sim
% time, so this is belt-and-braces: it mainly stops the block re-emitting
% the full N x 13 table at every continuous major step.
try
    set_param(BLK, 'SampleTime', char(opt.SampleTime));
    fprintf('  SampleTime set to %s s\n', char(opt.SampleTime));
catch err
    fprintf(['  NOTE: could not set SampleTime on the MATLAB Function ' ...
             'block (%s).\n        The 1 Hz gate inside ' ...
             'gfold_onboard_step still applies.\n'], err.identifier);
end

%% ---- wire it up ---------------------------------------------------------
fprintf('Wiring ...\n');

% inputs: t, x, m, xref_stored, unom_stored, tgrid
wire(mdl, sel.clock, 1, BLK, 1);
wire(mdl, sel.state, 1, BLK, 2);
wire(mdl, sel.mass,  1, BLK, 3);
wire(mdl, sel.xref,  1, BLK, 4);
wire(mdl, sel.unom,  1, BLK, 5);
wire(mdl, sel.tgrid, 1, BLK, 6);

% outputs: replanned tables into the Gain_Selector
wire(mdl, BLK, 1, sel.gain, 4);
wire(mdl, BLK, 2, sel.gain, 7);

% FDIR F-23 flag out to the base workspace
add_block('simulink/Sinks/To Workspace', FLAG, ...
    'Position', [gpos(1)-100, gpos(4)+40, gpos(1)-40, gpos(4)+70], ...
    'VariableName', 'gfold_flag', 'SaveFormat', 'Timeseries');
wire(mdl, BLK, 3, FLAG, 1);

%% ---- save --------------------------------------------------------------
if opt.Save
    save_system(mdl);
    fprintf('Saved %s\n', mdl_file);
else
    fprintf('Save skipped (Save=false); model left dirty in memory.\n');
end

fprintf(['\nDone. Run ctv/sim_setup_ctv.m to populate the reference and ' ...
         'gain tables,\nthen simulate %s.\n'], mdl);

end

% ========================================================================
function sel = locate(mdl)
% Find blocks by what they DO, not by their auto-generated names, so the
% script survives someone renaming Constant4..Constant8.

% Gain_Selector = the root subsystem that has an inport called K1flat
subs = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'SubSystem');
sel.gain = '';
for i = 1:numel(subs)
    inports = find_system(subs{i}, 'SearchDepth', 1, 'BlockType', 'Inport');
    nms = cellfun(@(b) get_param(b, 'Name'), inports, 'UniformOutput', false);
    if any(strcmp(nms, 'K1flat')) && any(strcmp(nms, 'xref'))
        sel.gain = subs{i};
        break
    end
end
assert(~isempty(sel.gain), ...
    'Could not find the Gain_Selector subsystem (an inport named K1flat).');

sel.tgrid = const_by_value(mdl, {'tgrid'});
% v2 names the reference table Target_Trajectory, v3 names it
% Target_Trajectory_full. Accept whichever this model actually uses.
sel.xref  = const_by_value(mdl, {'Target_Trajectory', 'Target_Trajectory_full'});
sel.unom  = const_by_value(mdl, {'unom'});

clk = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'Clock');
assert(~isempty(clk), 'Could not find a Clock block at the model root.');
sel.clock = clk{1};

% the 13-state vector is whatever currently feeds Gain_Selector port 2
sel.state = upstream(sel.gain, 2);

% mass = the root Saturate block clamped between wet and dry mass
sats = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'Saturate');
sel.mass = '';
for i = 1:numel(sats)
    ul = get_param(sats{i}, 'UpperLimit');
    if contains(ul, 'WetMass')
        sel.mass = sats{i};
        break
    end
end
assert(~isempty(sel.mass), ...
    ['Could not find the vehicle mass signal (a root Saturate block whose ' ...
     'UpperLimit mentions WetMass). Wire port 3 of GFOLD Replan by hand.']);

end

function blk = const_by_value(mdl, vals)
% vals is a cellstr of acceptable Value expressions, tried in order.
cs = find_system(mdl, 'SearchDepth', 1, 'BlockType', 'Constant');
for v = 1:numel(vals)
    for i = 1:numel(cs)
        if strcmp(strtrim(get_param(cs{i}, 'Value')), vals{v})
            blk = cs{i};
            return
        end
    end
end
error('Could not find a root Constant block with Value one of: %s.', ...
    strjoin(vals, ', '));
end

function src = upstream(blk, port)
ph = get_param(blk, 'PortHandles');
lh = get_param(ph.Inport(port), 'Line');
assert(lh ~= -1, 'Nothing is connected to port %d of %s.', port, blk);
sbh = get_param(lh, 'SrcBlockHandle');
assert(sbh ~= -1, 'Port %d of %s has a dangling line.', port, blk);
src = getfullname(sbh);
end

function chart = emchart(blkpath)
% Stateflow.EMChart lookup for a MATLAB Function block, with a fallback in
% case .Path is reported as the parent system rather than the block path.
rt = sfroot;
c  = rt.find('-isa', 'Stateflow.EMChart', 'Path', blkpath);
if isempty(c)
    [parent, name] = fileparts(blkpath);
    all_c = rt.find('-isa', 'Stateflow.EMChart');
    for i = 1:numel(all_c)
        if strcmp(all_c(i).Name, name) && strcmp(all_c(i).Path, parent)
            c = all_c(i);
            break
        end
    end
end
assert(~isempty(c), 'Could not get the Stateflow chart for %s', blkpath);
chart = c(1);
end

function drop_line(blk, port)
ph = get_param(blk, 'PortHandles');
lh = get_param(ph.Inport(port), 'Line');
if lh ~= -1
    delete_line(lh);
    fprintf('  removed existing line into %s port %d\n', ...
        get_param(blk, 'Name'), port);
end
end

function wire(mdl, srcBlk, srcPort, dstBlk, dstPort)
sp = get_param(srcBlk, 'PortHandles');
dp = get_param(dstBlk, 'PortHandles');
add_line(mdl, sp.Outport(srcPort), dp.Inport(dstPort), 'autorouting', 'on');
end

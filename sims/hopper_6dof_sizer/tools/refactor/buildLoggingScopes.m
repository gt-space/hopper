function buildLoggingScopes(model)
%BUILDLOGGINGSCOPES  Switchable scope groups and a full signal dump in the harness Logging block.
%   buildLoggingScopes()        models/hopper_6dof.slx (saves it)
%   buildLoggingScopes(model)   another loaded copy (saves it)
%
%   Logging gets every group's whole bus (env, nav, gnc, act, prop, dyn).
%   The loose scopes it had (X, Y, Z, Vz, Pitch, Scope6, CG Accel) give way
%   to one Variant Subsystem per group in scopeGroups, Scopes_<name>,
%   active when VARIANT_SCOPES_<NAME> == 2, and Signal_Dump, active when
%   VARIANT_LOG == 2, which logs all six buses (every signal between the
%   groups) to the simulation's logsout. variantConfig('SCOPES', ...,
%   'LOG', ...) sets them; both default to off, and choices that are off
%   aren't compiled, so the default flight is unchanged. The To Workspace
%   blocks the tests and Monte Carlo read stay as they are.
%
%   Rerunning rebuilds the groups from scopeGroups, so editing that list
%   and rerunning is how a group changes.

if nargin < 1, model = 'hopper_6dof'; end
load_system(model);
L = [model '/Logging'];
buses = {'env_bus', 'Environment'; 'nav_bus', 'Navigation_Sensors'; 'gnc_bus', 'GNC'
         'act_bus', 'Actuation'; 'prop_bus', 'Mass_Propulsion'; 'dyn_bus', 'Vehicle_Dynamics'};

% 1. the whole of each bus into Logging (an In Bus Element with no element)
for k = 1:size(buses, 1)
    all_ = [buses{k, 1} '_all'];   % a block can't take its port's name
    if getSimulinkBlockHandle([L '/' all_]) > 0, continue, end
    old = find_system(L, 'SearchDepth', 1, 'BlockType', 'Inport', 'IsBusElementPort', 'on', 'PortName', buses{k, 1});
    if isempty(old)
        b = add_block('simulink/Ports & Subsystems/In Bus Element', [L '/' all_]);
        set_param(b, 'PortName', buses{k, 1});
        set_param(b, 'Element', '');
        sp = get_param([model '/' buses{k, 2}], 'PortHandles');
        dp = get_param(L, 'PortHandles');
        add_line(model, sp.Outport(1), dp.Inport(str2double(get_param(b, 'Port'))), 'autorouting', 'smart');
    else
        b = add_block(old{1}, [L '/' all_], 'CopyOption', 'duplicate');
        set_param(b, 'Element', '');
    end
end

% 2. the old loose scopes and what fed only them
deleteBlocks(L, [strcat([L '/'], {'CG Accel', 'Pitch', 'Scope6', 'Vz', 'X', 'Y', 'Z', ...
    'From Workspace', 'From Workspace1', 'From Workspace2', 'Radians to Degrees1'}), ...
    find_system(L, 'SearchDepth', 1, 'RegExp', 'on', 'Name', '^(Scopes_|Signal_Dump$)')']);
for b = find_system(L, 'SearchDepth', 1, 'BlockType', 'Inport', 'IsBusElementPort', 'on')'
    p = get_param(b{1}, 'PortHandles');
    if ~isempty(get_param(b{1}, 'Element')) && get_param(p.Outport, 'Line') < 0, delete_block(b{1}); end
end

% 3. one switchable scope group per scopeGroups entry, and the dump
G = scopeGroups();
for g = 1:numel(G)
    makeSwitch(L, ['Scopes_' G(g).name], G(g).buses, sprintf('VARIANT_SCOPES_%s == 2', upper(G(g).name)), ...
        @(c) scopeChoice(c, G(g)));
end
makeSwitch(L, 'Signal_Dump', buses(:, 1)', 'VARIANT_LOG == 2', @(c) dumpChoice(c, buses(:, 1)'));

layoutLogging(L, {G.name});
save_system(model);
fprintf('%s: %d scope groups and the signal dump in Logging\n', model, numel(G));
end

function layoutLogging(L, names)
% Two columns: on the left one row per To Workspace log (its bus element,
% any conversion, the block); on the right the whole buses and the scope
% groups and dump they feed
y = 30;
tw = find_system(L, 'SearchDepth', 1, 'BlockType', 'ToWorkspace');
[~, o] = sort(lower(get_param(tw, 'VariableName'))); tw = tw(o);
for k = 1:numel(tw)
    chain = tw(k);   % back from the log to its bus element
    while true
        p = get_param(chain{1}, 'PortHandles');
        if isempty(p.Inport) || get_param(p.Inport(1), 'Line') < 0, break, end
        chain = [{getfullname(get_param(get_param(p.Inport(1), 'Line'), 'SrcBlockHandle'))}, chain]; %#ok<AGROW>
    end
    xs = [40, 200, 300]; w = [20, 50, 160]; hs = [14, 30, 30];   % element, conversion, log
    col = [1, 3]; if numel(chain) == 3, col = [1 2 3]; end
    for j = 1:numel(chain)
        c = col(j);
        set_param(chain{j}, 'Orientation', 'right');   % some were flipped or rotated
        if c == 3, set_param(chain{j}, 'ShowName', 'off'); end   % the log's variable says it
        yc = y + (30 - hs(c)) / 2;   % centred on the row, so lines run straight
        set_param(chain{j}, 'Position', [xs(c), yc, xs(c) + w(c), yc + hs(c)]);
    end
    y = y + 50;
end
x0 = 620; y = 30;
for n = [strcat('Scopes_', names), {'Signal_Dump'}]
    b = [L '/' n{1}]; np = get_param(b, 'Ports');
    h = 30 * np(1) + 20;
    set_param(b, 'Position', [x0 + 260, y, x0 + 420, y + h]);
    y = y + h + 40;
end
ins = find_system(L, 'SearchDepth', 1, 'RegExp', 'on', 'Name', '_bus_all$');
for k = 1:numel(ins)
    yy = 30 + (y - 60) * (k - 0.5) / numel(ins);
    set_param(ins{k}, 'Position', [x0, yy, x0 + 20, yy + 14]);
end
Simulink.BlockDiagram.routeLine(find_system(L, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line'));
end

function makeSwitch(L, name, buses, cond, build)
% Variant Subsystem L/name fed by the whole buses, with one choice, On,
% active when cond holds and otherwise nothing at all
blk = [L '/' name];
add_block('simulink/Ports & Subsystems/Variant Subsystem', blk);
Simulink.SubSystem.deleteContents(blk);
for k = 1:numel(buses)
    add_block('built-in/Inport', [blk '/' buses{k}], 'Position', [20, 40 * k, 50, 40 * k + 14]);
end
c = [blk '/On'];
add_block('built-in/Subsystem', c, 'Position', [200 60 360 160]);
build(c);
set_param(c, 'VariantControl', cond);
set_param(blk, 'VariantControlMode', 'expression', 'VariantActivationTime', 'update diagram', ...
    'PropagateVariantConditions', 'off', 'AllowZeroVariantControls', 'on');
bp = get_param(blk, 'PortHandles');
for k = 1:numel(buses)
    sp = get_param([L '/' buses{k} '_all'], 'PortHandles');
    add_line(L, sp.Outport, bp.Inport(k), 'autorouting', 'smart');
end
end

function scopeChoice(c, g)
% Inports -> one Bus Selector per bus -> transforms -> one multi-display scope
sel = containers.Map();
sig = [g.displays{:, 2}];
for k = 1:numel(g.buses)
    b = g.buses{k};
    add_block('built-in/Inport', [c '/' b]);
    el = unique(cellfun(@(s) extractBefore([extractAfter(s, '.') ':'], ':'), ...
        sig(startsWith(sig, [b '.'])), 'uni', 0), 'stable');
    bs = add_block('simulink/Signal Routing/Bus Selector', [c '/sel_' b], 'OutputSignals', strjoin(el, ','));
    add_line(c, [b '/1'], ['sel_' b '/1']);
    op = get_param(bs, 'PortHandles');
    for e = 1:numel(el), sel([b '.' el{e}]) = op.Outport(e); end
end
n = size(g.displays, 1);
sc = add_block('simulink/Sinks/Scope', [c '/' g.name]);
cfg = get_param(sc, 'ScopeConfiguration');
cfg.NumInputPorts = num2str(n);
cfg.LayoutDimensions = [n 1];
cfg.OpenAtSimulationStart = true;
for d = 1:n
    cfg.ActiveDisplay = d;
    cfg.Title = g.displays{d, 1};
    cfg.ShowLegend = true;
end
scp = get_param(sc, 'PortHandles');
for d = 1:n
    s = g.displays{d, 2}; src = zeros(1, numel(s)); nm = cell(1, numel(s));
    for k = 1:numel(s), [src(k), nm{k}] = signal(c, s{k}, sel); end
    if numel(src) == 1
        nameLine(add_line(c, src, scp.Inport(d), 'autorouting', 'smart'), nm{1});
    else
        mx = add_block('simulink/Signal Routing/Mux', [c '/mux'], 'MakeNameUnique', 'on', 'Inputs', num2str(numel(src)));
        mp = get_param(mx, 'PortHandles');
        for k = 1:numel(src), nameLine(add_line(c, src(k), mp.Inport(k), 'autorouting', 'smart'), nm{k}); end
        add_line(c, mp.Outport, scp.Inport(d), 'autorouting', 'smart');
    end
end
try Simulink.BlockDiagram.arrangeSystem(c); catch, end
end

function [p, name] = signal(c, s, sel)
% the output port carrying signal spec s (see scopeGroups), and a name for its line
if startsWith(s, 'ref.')
    v = extractAfter(s, '.');
    fw = add_block('simulink/Sources/From Workspace', [c '/' v], 'MakeNameUnique', 'on', 'VariableName', v);
    p = get_param(fw, 'PortHandles').Outport;
    name = v;
else
    [base, tf] = strtok(s, ':'); tf = erase(tf, ':');
    p = sel(base);
    name = extractAfter(base, '.');
    switch extractBefore([tf '0'], 2)   % first letter: d(eg/iag), i(ndices)
        case 'd'
            if strcmp(tf, 'deg')
                p = chain(c, p, 'simulink/Math Operations/Gain', {'Gain', '180/pi'});
                name = [name '_deg'];
            else   % diag
                p = chain(c, p, 'simulink/Math Operations/Reshape', {'OutputDimensionality', '1-D array'});
                p = chain(c, p, 'simulink/Signal Routing/Selector', {'Indices', '[1 5 9]', 'InputPortWidth', '9'});
                name = [name '_diag'];
            end
        case 'i'
            r = sscanf(tf, 'i%d-%d');
            p = chain(c, p, 'simulink/Signal Routing/Selector', {'Indices', sprintf('[%d:%d]', r(1), r(2)), 'InputPortWidth', '-1'});
            name = sprintf('%s_%d_%d', name, r(1), r(2));
    end
end
end

function nameLine(l, name)
% label a scope input for the legend; a Bus Selector's outputs already
% carry their element names and can't be relabelled
if ~strcmp(get_param(get_param(l, 'SrcBlockHandle'), 'BlockType'), 'BusSelector'), set_param(l, 'Name', name); end
end

function p = chain(c, src, lib, params)
b = add_block(lib, [c '/' regexprep(lib, '.*/', '')], 'MakeNameUnique', 'on', params{:});
bp = get_param(b, 'PortHandles');
add_line(c, src, bp.Inport(1), 'autorouting', 'smart');
p = bp.Outport(1);
end

function dumpChoice(c, buses)
% each bus logged whole to logsout under its own name, then terminated
for k = 1:numel(buses)
    in = add_block('built-in/Inport', [c '/' buses{k}], 'Position', [20, 50 * k, 50, 50 * k + 14]);
    add_block('built-in/Terminator', [c '/end_' buses{k}], 'Position', [200, 50 * k, 220, 50 * k + 14]);
    add_line(c, [buses{k} '/1'], ['end_' buses{k} '/1']);
    p = get_param(in, 'PortHandles');
    set_param(p.Outport, 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', buses{k});
end
end

function deleteBlocks(sys, doomed)
% Delete blocks in sys, keeping every connection among the others (a
% delete_line on a branched line takes the whole tree, so redraw it)
doomed = doomed(cellfun(@(b) getSimulinkBlockHandle(b) > 0, doomed));
if isempty(doomed), return, end
dh = cellfun(@getSimulinkBlockHandle, doomed);
srcs = [];   % source ports whose line trees touch a doomed block
for h = dh
    p = get_param(h, 'PortHandles');
    srcs = [srcs, p.Outport]; %#ok<AGROW>
    for i = p.Inport
        l = get_param(i, 'Line'); if l > 0, srcs(end+1) = get_param(l, 'SrcPortHandle'); end %#ok<AGROW>
    end
end
srcs = unique(srcs(srcs > 0));
keep = cell(0, 2);
for b = find_system(sys, 'SearchDepth', 1, 'Type', 'Block')'
    h = getSimulinkBlockHandle(b{1});
    if strcmp(b{1}, sys) || any(h == dh), continue, end
    p = get_param(h, 'PortHandles');
    for i = p.Inport
        l = get_param(i, 'Line');
        if l > 0
            s = get_param(l, 'SrcPortHandle');
            if any(s == srcs) && ~any(get_param(s, 'ParentHandle') == dh), keep(end+1, :) = {s, i}; end %#ok<AGROW>
        end
    end
end
for s = srcs
    l = get_param(s, 'Line'); if l > 0, delete_line(l); end
end
for h = dh, delete_block(h); end
for k = 1:size(keep, 1), add_line(sys, keep{k, 1}, keep{k, 2}, 'autorouting', 'smart'); end
end

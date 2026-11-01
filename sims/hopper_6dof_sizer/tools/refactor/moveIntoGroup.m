function grp = moveIntoGroup(model, group, members)
%MOVEINTOGROUP  Move top-level blocks into a new virtual subsystem, rewired explicitly.
%   grp = moveIntoGroup(model, group, members)
%     model    loaded model
%     group    name of the new subsystem at the top level
%     members  top-level block names (newlines written as spaces)
%
%   Unlike Simulink.BlockDiagram.createSubsystem (which changed the
%   harness's flight and crashed MATLAB), every connection is recorded by
%   block name and port number first, then rebuilt:
%     - inside the group: member to member, redrawn inside
%     - coming in: the outside source gets a new branch to a group inport
%       (one inport per distinct source port); live lines are never deleted
%     - going out: one group outport per member output, wired to every
%       outside destination; only these lines, whose source block is gone,
%       are deleted
%   Ports are named after the signal they carry. The subsystem is virtual,
%   so the simulation is unchanged.

grp = [model '/' group];
assert(getSimulinkBlockHandle(grp) < 0, 'moveIntoGroup:exists', '%s already exists', grp);
path = cellfun(@(n) topBlock(model, n), members, 'UniformOutput', false);
isMember = @(blk) any(strcmp(path, blk));

% 1. record every connection of the members
inc = cell(0, 4); internal = cell(0, 4); out = cell(0, 4);   % rows: {srcBlk, srcPort, dstBlk, dstPort}
for k = 1:numel(path)
    ph = get_param(path{k}, 'PortHandles');
    for i = 1:numel(ph.Inport)
        l = get_param(ph.Inport(i), 'Line');
        if l < 0, continue, end
        s = get_param(l, 'SrcPortHandle');
        if s < 0, continue, end
        sb = get_param(s, 'Parent'); sp = get_param(s, 'PortNumber');
        if isMember(sb)
            internal(end+1, :) = {name(sb), sp, name(path{k}), i}; %#ok<AGROW>
        else
            inc(end+1, :) = {sb, sp, name(path{k}), i}; %#ok<AGROW>
        end
    end
    for o = 1:numel(ph.Outport)
        l = get_param(ph.Outport(o), 'Line');
        if l < 0, continue, end
        d = get_param(l, 'DstPortHandle');
        for dh = d(d > 0)'
            db = get_param(dh, 'Parent');
            if ~isMember(db)
                out(end+1, :) = {name(path{k}), o, db, get_param(dh, 'PortNumber')}; %#ok<AGROW>
            end
        end
    end
end
for k = 1:numel(path)
    ph = get_param(path{k}, 'PortHandles');
    assert(isempty([ph.Enable, ph.Trigger, ph.LConn, ph.RConn]), 'moveIntoGroup:ports', ...
        '%s has enable, trigger or physical ports, which this does not move', path{k});
end

% 2. an empty subsystem; copy the members in, keeping their relative layout
pos = cell2mat(get_param(path, 'Position'));
add_block('built-in/Subsystem', grp, 'Position', [min(pos(:, 1)), min(pos(:, 2)), min(pos(:, 1)) + 200, min(pos(:, 2)) + 200]);
for k = 1:numel(path)
    add_block(path{k}, [grp '/' name(path{k})], 'Position', get_param(path{k}, 'Position') - [pos(1, 1:2) pos(1, 1:2)] + 300);
end

% 3. outgoing lines lose their source with the originals: delete them
for k = 1:numel(path)
    ph = get_param(path{k}, 'PortHandles');
    for o = ph.Outport
        l = get_param(o, 'Line');
        if l > 0, delete_line(l); end
    end
    delete_block(path{k});
end

% 4. inside
for r = 1:size(internal, 1)
    add_line(grp, sprintf('%s/%d', internal{r, 1}, internal{r, 2}), sprintf('%s/%d', internal{r, 3}, internal{r, 4}), 'autorouting', 'smart');
end

% 5. coming in: one inport per distinct outside source port
names = {};
[keys, ~, which] = unique(cellfun(@(b, p) sprintf('%s|%d', b, p), inc(:, 1), inc(:, 2), 'UniformOutput', false), 'stable');
for k = 1:numel(keys)
    rows = find(which == k);
    sb = inc{rows(1), 1}; sp = inc{rows(1), 2};
    [n, names] = freeName(signalName(sb, sp), names, members);
    add_block('built-in/Inport', [grp '/' n], 'Position', [20, 40 * k, 50, 40 * k + 14]);
    for r = rows'
        add_line(grp, [n '/1'], sprintf('%s/%d', inc{r, 3}, inc{r, 4}), 'autorouting', 'smart');
    end
    gp = get_param(grp, 'PortHandles');
    src = get_param(sb, 'PortHandles');
    l = add_line(model, src.Outport(sp), gp.Inport(k), 'autorouting', 'smart');
    set_param(l, 'Name', n);
end

% 6. going out: one outport per member output port
[keys, ~, which] = unique(cellfun(@(b, p) sprintf('%s|%d', b, p), out(:, 1), out(:, 2), 'UniformOutput', false), 'stable');
for k = 1:numel(keys)
    rows = find(which == k);
    mb = out{rows(1), 1}; mp = out{rows(1), 2};
    [n, names] = freeName(signalName([grp '/' mb], mp), names, members);
    add_block('built-in/Outport', [grp '/' n], 'Position', [900, 40 * k, 930, 40 * k + 14]);
    add_line(grp, sprintf('%s/%d', mb, mp), [n '/1'], 'autorouting', 'smart');
    gp = get_param(grp, 'PortHandles');
    for r = rows'
        dst = get_param(out{r, 3}, 'PortHandles');
        add_line(model, gp.Outport(k), dst.Inport(out{r, 4}), 'autorouting', 'smart');
    end
end
gp = get_param(grp, 'PortHandles');
fprintf('%s: %d blocks, %d in, %d out, %d inside\n', grp, numel(path), numel(gp.Inport), numel(gp.Outport), size(internal, 1));
end

function p = topBlock(model, n)
b = find_system(model, 'SearchDepth', 1, 'Type', 'Block');
k = find(strcmp(strrep(b, newline, ' '), [model '/' n]));
assert(isscalar(k), 'moveIntoGroup:block', 'no unique top-level block %s', n);
p = b{k};
end

function n = name(blk)
n = get_param(blk, 'Name');
end

function n = signalName(blk, port)
% the source port's own name if it has one, else block (and port number)
n = '';
if strcmp(get_param(blk, 'BlockType'), 'SubSystem')
    o = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Outport', 'Port', num2str(port));
    if ~isempty(o), n = get_param(o{1}, 'Name'); end
end
b = clean(get_param(blk, 'Name'));
p = get_param(blk, 'Ports');
if isempty(n) || ~isempty(regexp(n, '^Out\d*$', 'once'))
    n = b;
    if p(2) > 1, n = sprintf('%s_%d', b, port); end
else
    n = [b '_' clean(n)];
end
n = clean(n);
end

function s = clean(s)
s = regexprep(strtrim(strrep(s, newline, ' ')), '[^\w]+', '_');
s = regexprep(s, '^_+|_+$', '');
end

function [n, names] = freeName(n, names, members)
taken = [names, cellfun(@clean, members, 'UniformOutput', false), members];
base = n; k = 1;
while any(strcmp(taken, n)), k = k + 1; n = sprintf('%s_%d', base, k); end
names{end+1} = n;
end

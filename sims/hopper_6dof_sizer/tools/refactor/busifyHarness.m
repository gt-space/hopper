function busifyHarness(model)
%BUSIFYHARNESS  Connect the harness groups with one bus per group instead of ~100 wires.
%   busifyHarness()        models/hopper_6dof.slx (saves it)
%   busifyHarness(model)   another loaded copy (saves it)
%
%   After moveIntoGroup, every top-level connection runs group port to
%   group port. Here each group's outports become Out Bus Element blocks
%   on one bus port (env_bus, nav_bus, ...); each group inport becomes an
%   In Bus Element block picking its element from the source group's bus.
%   Every element gets a short snake_case name from signalNames() below,
%   and the bus element blocks carry that name too. The top level is then
%   one line per source-destination pair of groups. Virtual buses don't
%   change the simulation. Rerunning is harmless.

if nargin < 1, model = 'hopper_6dof'; end
load_system(model);
if ~isempty(find_system(model, 'SearchDepth', 2, 'BlockType', 'Outport', 'IsBusElementPort', 'on'))
    fprintf('%s: already on buses\n', model); return
end
[busOf, names] = signalNames();
grp = find_system(model, 'SearchDepth', 1, 'BlockType', 'SubSystem');
grp = grp(ismember(get_param(grp, 'Name'), [keys(busOf) {'Logging'}]));
deleteStrays(model, grp);

% 1. element name of each group output, from what drives the inside Outport
elem = containers.Map();
for g = 1:numel(grp)
    G = grp{g};
    for h = find_system(G, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'Outport')'
        p = get_param(h, 'PortHandles');
        s = get_param(get_param(p.Inport, 'Line'), 'SrcPortHandle');
        key = sprintf('%s|%s:%d', blockName(G), clean(get_param(get_param(s, 'Parent'), 'Name')), get_param(s, 'PortNumber'));
        assert(isKey(names, key), 'busifyHarness:unnamed', 'no signal name for %s', key);
        elem(sprintf('%s|%s', G, get_param(h, 'Port'))) = names(key);
    end
end
% what feeds each group inport: {dstGroup, dstPort, srcGroup, srcPort}
conn = cell(0, 4);
for d = 1:numel(grp)
    ph = get_param(grp{d}, 'PortHandles');
    for j = 1:numel(ph.Inport)
        l = get_param(ph.Inport(j), 'Line');
        assert(l > 0, 'busifyHarness:open', '%s input %d is unconnected', grp{d}, j);
        s = get_param(l, 'SrcPortHandle');
        sb = get_param(s, 'Parent');
        assert(any(strcmp(grp, sb)), 'busifyHarness:source', '%s input %d is fed by %s, not a group', grp{d}, j, sb);
        conn(end+1, :) = {grp{d}, j, sb, get_param(s, 'PortNumber')}; %#ok<AGROW>
    end
end

% 2. inside each group: outports -> Out Bus Element on the group's bus
for g = 1:numel(grp)
    G = grp{g}; first = [];
    o = find_system(G, 'SearchDepth', 1, 'FindAll', 'on', 'BlockType', 'Outport');
    on = arrayfun(@(h) elem(sprintf('%s|%s', G, get_param(h, 'Port'))), o, 'UniformOutput', false);   % before deletes renumber ports
    for k = 1:numel(o)
        h = o(k); n = on{k};
        pos = get_param(h, 'Position'); p = get_param(h, 'PortHandles');
        src = get_param(get_param(p.Inport, 'Line'), 'SrcPortHandle');
        delete_line(G, portRef(src), [esc(get_param(h, 'Name')) '/1']);   % only this branch: the source may feed others
        delete_block(h);
        [first, b] = addBusElement(G, 'Out', first, busOf(blockName(G)), n, pos);
        bp = get_param(b, 'PortHandles');
        add_line(G, src, bp.Inport, 'autorouting', 'smart');
    end
end

% 3. inside each group: inports -> In Bus Element from the source's bus
for r = 1:size(conn, 1)
    D = conn{r, 1}; j = conn{r, 2};
    in = find_system(D, 'SearchDepth', 1, 'BlockType', 'Inport', 'Port', num2str(j));
    p = get_param(in{1}, 'PortHandles');
    l = get_param(p.Outport, 'Line');
    dst = []; if l > 0, dst = get_param(l, 'DstPortHandle'); dst = dst(dst > 0); delete_line(l); end   % all its own branches
    conn{r, 5} = in{1}; conn{r, 6} = elem(sprintf('%s|%d', conn{r, 3}, conn{r, 4})); %#ok<AGROW>
    conn{r, 7} = dst; conn{r, 8} = get_param(in{1}, 'Position'); %#ok<AGROW>
end
% delete the old inports only after every one is recorded: port numbers shift on delete
for r = 1:size(conn, 1), delete_block(conn{r, 5}); end
firstIn = containers.Map();   % '<dst group>|<bus>' -> first In Bus Element on that port
for r = 1:size(conn, 1)
    D = conn{r, 1}; bus = busOf(blockName(conn{r, 3})); k = [D '|' bus];
    f = []; if isKey(firstIn, k), f = firstIn(k); end
    [f, b] = addBusElement(D, 'In', f, bus, conn{r, 6}, conn{r, 8}); firstIn(k) = f;
    bp = get_param(b, 'PortHandles');
    for d = conn{r, 7}', add_line(D, bp.Outport, d, 'autorouting', 'smart'); end
end

% 4. the top level: every old line goes (stubs included); one bus line per pair
delete_line(find_system(model, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'line'));
pairs = unique(cellfun(@(s, d) [s '|' d], conn(:, 3), conn(:, 1), 'UniformOutput', false));
for k = 1:numel(pairs)
    sd = strsplit(pairs{k}, '|'); S = sd{1}; D = sd{2}; bus = busOf(blockName(S));
    sp = get_param(S, 'PortHandles');
    ibe = find_system(D, 'SearchDepth', 1, 'BlockType', 'Inport', 'IsBusElementPort', 'on', 'PortName', bus);
    dp = get_param(D, 'PortHandles');
    add_line(model, sp.Outport(1), dp.Inport(str2double(get_param(ibe{1}, 'Port'))), 'autorouting', 'smart');
end
% each bus is labelled once, on the line leaving its group (branches share it)
for g = 1:numel(grp)
    sp = get_param(grp{g}, 'PortHandles');
    if ~isempty(sp.Outport), set_param(get_param(sp.Outport(1), 'Line'), 'Name', busOf(blockName(grp{g}))); end
end
% the top level is only groups and buses now: notes left there go
delete(find_system(model, 'SearchDepth', 1, 'FindAll', 'on', 'Type', 'annotation'));
save_system(model);
fprintf('%s: %d group connections now %d bus lines\n', model, size(conn, 1), numel(pairs));
end

function [busOf, names] = signalNames()
% Bus per group, and the element name of every group output, keyed by
% '<group>|<inside source block>:<port>'. Short, snake_case, unique.
busOf = containers.Map( ...
    {'Environment', 'GNC', 'Actuation', 'Mass_Propulsion', 'Navigation_Sensors', 'Vehicle_Dynamics'}, ...
    {'env_bus', 'gnc_bus', 'act_bus', 'prop_bus', 'nav_bus', 'dyn_bus'});
t = {
    'Environment|MatrixMultiply1:1'                 'wind_force'
    'Environment|MatrixMultiply2:1'                 'wind_moment'
    'Environment|Product1:1'                        'roll_dist'
    'GNC|Actuator_Faults:1'                         'u_cmd'
    'GNC|FDIR:2'                                    'nav_mask'
    'GNC|FDIR:3'                                    'nav_r_scale'
    'Actuation|Actuators_TVC_RCS:1'                 'force'
    'Actuation|Actuators_TVC_RCS:2'                 'moment'
    'Actuation|Actuators_TVC_RCS:3'                 'tvc_pitch'
    'Actuation|Actuators_TVC_RCS:4'                 'tvc_yaw'
    'Actuation|Quantizer:1'                         'thrust_cmd'
    'Actuation|IC:1'                                'thrust_req'
    'Mass_Propulsion|Propulsion:1'                  'thrust'
    'Mass_Propulsion|Propellant_Mass_Integration:1' 'ox_mass'
    'Mass_Propulsion|Propellant_Mass_Integration:2' 'fu_mass'
    'Mass_Propulsion|Saturation1:1'                 'mass'
    'Mass_Propulsion|Mass_Properties:1'             'cg'
    'Mass_Propulsion|Mass_Properties:2'             'moi'
    'Mass_Propulsion|Mass_Properties:3'             'moi_dot'
    'Mass_Propulsion|Mass_Properties:4'             'engine_cg'
    'Mass_Propulsion|Mass_Properties:5'             'cg_slosh'
    'Mass_Propulsion|Mass_Properties:6'             'moi_slosh'
    'Navigation_Sensors|Sensors:1'                  'meas'
    'Navigation_Sensors|Sensors:2'                  'valid'
    'Navigation_Sensors|Sensors:3'                  'stamp'
    'Navigation_Sensors|Navigation:1'               'x_nav'
    'Navigation_Sensors|Navigation:2'               'dcm_nav'
    'Navigation_Sensors|Navigation:3'               'x_filt'
    'Navigation_Sensors|Navigation:4'               'p_cov'
    'Navigation_Sensors|Navigation:5'               'nav_diag'
    'Navigation_Sensors|Navigation:6'               'shadow'
    'Navigation_Sensors|Navigation:7'               'x_vert'
    'Vehicle_Dynamics|6DOF (Quaternion):1'          'vel_ned'
    'Vehicle_Dynamics|6DOF (Quaternion):4'          'quat'
    'Vehicle_Dynamics|6DOF (Quaternion):5'          'dcm_be'
    'Vehicle_Dynamics|6DOF (Quaternion):7'          'w_body'
    'Vehicle_Dynamics|Mux2:1'                       'pos_ned'
    'Vehicle_Dynamics|Sum11:1'                      'x'
    'Vehicle_Dynamics|Sum14:1'                      'y'
    'Vehicle_Dynamics|Sum15:1'                      'z'
    'Vehicle_Dynamics|Demux1:1'                     'phi'
    'Vehicle_Dynamics|Demux1:2'                     'theta'
    'Vehicle_Dynamics|Demux1:3'                     'psi'
    'Vehicle_Dynamics|Demux3:1'                     'u'
    'Vehicle_Dynamics|Demux3:3'                     'w'
    'Vehicle_Dynamics|Demux4:1'                     'p'
    'Vehicle_Dynamics|Slosh Specific Force:1'       'f_body'
    'Vehicle_Dynamics|Sum6:1'                       'cg_accel'
    'Vehicle_Dynamics|MatrixMultiply:1'             'g_body'
    'Vehicle_Dynamics|Slosh:3'                      'ox_slosh_pos'
    'Vehicle_Dynamics|Slosh:4'                      'ox_slosh_mass'
    'Vehicle_Dynamics|Slosh:5'                      'fu_slosh_pos'
    'Vehicle_Dynamics|Slosh:6'                      'fu_slosh_mass'
    };
assert(numel(unique(t(:, 2))) == size(t, 1), 'busifyHarness:dupName', 'signal names must be unique');
names = containers.Map(t(:, 1), t(:, 2));
end

function [first, b] = addBusElement(G, kind, first, bus, n, pos)
% One bus element block named n on the bus port. The first on a port
% comes from the library with its PortName; later ones duplicate it onto
% the same port (a new block can't take a port name that exists).
name = freeName(G, n, lower(kind));
if isempty(first)
    b = add_block(['simulink/Ports & Subsystems/' kind ' Bus Element'], [G '/' esc(name)], 'Position', pos);
    set_param(b, 'PortName', bus, 'Element', n);
    first = b;
else
    b = add_block(getfullname(first), [G '/' esc(name)], 'CopyOption', 'duplicate', 'Position', pos);
    set_param(b, 'Element', n);
end
end

function n = freeName(G, n, kind)
% n, or n_<kind> if some other block in the group already has that name
if getSimulinkBlockHandle([G '/' esc(n)]) > 0, n = [n '_' kind]; end
assert(getSimulinkBlockHandle([G '/' esc(n)]) < 0, 'busifyHarness:name', '%s/%s exists', G, n);
end

function deleteStrays(model, grp)
% top-level blocks outside the groups (left-over taps whose scopes are
% gone) must feed nothing; they would be cut loose by step 4 anyway
top = find_system(model, 'SearchDepth', 1, 'Type', 'Block');
for b = top(~ismember(top, [{model}; grp]))'
    ph = get_param(b{1}, 'PortHandles');
    for o = ph.Outport
        l = get_param(o, 'Line');
        assert(l < 0 || all(get_param(l, 'DstPortHandle') < 0), 'busifyHarness:stray', ...
            '%s still feeds something', b{1});
    end
    fprintf('deleting stray top-level block "%s"\n', clean(get_param(b{1}, 'Name')));
    delete_block(b{1});
end
end

function n = blockName(p)
n = get_param(p, 'Name');
end

function s = clean(s)
s = strrep(s, newline, ' ');
end

function s = esc(n)
s = strrep(n, '/', '//');
end

function r = portRef(ph)
r = sprintf('%s/%d', esc(get_param(get_param(ph, 'Parent'), 'Name')), get_param(ph, 'PortNumber'));
end

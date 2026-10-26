function build_hopper_6dof()
%BUILD_HOPPER_6DOF  Build the variant harness models/hopper_6dof.slx from v2.
%   Copies hopper_6dof_NED_v2 and puts Variant Subsystems at the points
%   where the feature copies of v2 differ, one per axis in variantAxes:
%
%     Slosh            SLOSH       full: subsystems/slosh/Slosh.slx
%                                  off:  zero slosh force, moment, offsets
%     Propulsion       PROPULSION  lut: thrust command passes through,
%                                       Thrust_to_Mdot_LUT gives the mdots
%                                  fluid_network: #26's Fluid_Network
%                                       (thrust and mdots from the feed system)
%     Propellant_Flow  MASS_PROPS  varying: mdots pass through
%                                  constant: zero mdot, so propellant mass,
%                                       CG and MOI hold their initial values
%                                       and dI/dt = 0
%     Ground_Contact   GROUND      none: zero force and moment (as v2)
%                                  strut: subsystems/ground/Ground_Contact.slx
%                                       into Sum4/Sum5, v2's unused inputs
%
%   Choice k of axis X is active when VARIANT_X == k (variantConfig sets
%   them). With every axis at its default the harness is v2. Also removes
%   v2's commented-out leftovers (the old Fluid_Network, the scope-only
%   Ground_Contact_Unwired). Run with the HopperSim project open, after
%   convertToSubsystemReferences and build_v2_fluid_integration.
%   Writes models/hopper_6dof.slx.

base = 'hopper_6dof_NED_v2';
fluid = 'hopper_6dof_NED_v2_fluid_integration';
m = 'hopper_6dof';
out = hopperRoot('models', [m '.slx']);

if bdIsLoaded(m), close_system(m, 0); end
load_system(base);
save_system(base, out);
close_system(base, 0);
load_system(fluid);
cleanup = onCleanup(@() close_system(fluid, 0));

% v2 leftovers that the variants replace
for b = {'Fluid_Network', 'Ground_Contact_Unwired'}
    deleteWithLines([m '/' b{1}]);
end

% SLOSH ------------------------------------------------------------------
sloshPorts = portNames([m '/Slosh']);
makeVariant(m, 'Slosh', 'SLOSH', {
    'Slosh_Full', @(c) add_block([m '/Slosh'], c)      % the existing reference
    'Slosh_Off',  @(c) zeroChoice(c, sloshPorts.in, sloshPorts.out, ...
        {'[0 0 0]', 'zeros(3,1)', '[0 0 0]', '[0 0 0]', '[0 0 0]', '[0 0 0]'})
    });

% PROPULSION -------------------------------------------------------------
% v2: Propellant_Depletion_Cutoff -> (actuators, logging, Thrust_to_Mdot_LUT)
lut = [m '/Thrust_to_Mdot_LUT'];
cutoff = [m '/Propellant_Depletion_Cutoff'];
thrustDst = destinations(cutoff, 1);
thrustDst = thrustDst(~strcmp(thrustDst(:, 1), 'Thrust_to_Mdot_LUT'), :);
mdotDst = {destinations(lut, 1), destinations(lut, 2)};   % LUT out 1 feeds the ox side
lutPos = get_param(lut, 'Position');
deleteWithLines(lut);
for k = 1:size(thrustDst, 1)
    disconnectInput([m '/' thrustDst{k, 1}], thrustDst{k, 2});
end
disconnectInput([m '/IC1'], 1);
props = {{'Thrust_Cmd'}; {'Thrust'; 'Ox_mdot'; 'Fu_mdot'}};
add_block('built-in/Subsystem', [m '/Propulsion'], 'Position', lutPos);
makeVariant(m, 'Propulsion', 'PROPULSION', {
    'Propulsion_LUT', @(c) lutChoice(c, props)
    'Propulsion_Fluid_Network', @(c) fluidChoice(c, props, [fluid '/Fluid_Network'])
    });
connect(m, 'Propellant_Depletion_Cutoff', 1, 'Propulsion', 1);
connect(m, 'Propulsion', 1, 'IC1', 1);
for k = 1:size(thrustDst, 1)
    connect(m, 'Propulsion', 1, thrustDst{k, 1}, thrustDst{k, 2});
end

% MASS_PROPS: between the mdots and everything that reads them ----------
add_block('built-in/Subsystem', [m '/Propellant_Flow'], 'Position', lutPos + [120 0 120 0]);
flow = {{'Ox_mdot_in'; 'Fu_mdot_in'}; {'Ox_mdot'; 'Fu_mdot'}};
makeVariant(m, 'Propellant_Flow', 'MASS_PROPS', {
    'Flow_Varying', @(c) passChoice(c, flow)
    'Flow_Constant', @(c) zeroChoice(c, flow{1}, flow{2}, {'0', '0'})
    });
for o = 1:2
    connect(m, 'Propulsion', o + 1, 'Propellant_Flow', o);
    for k = 1:size(mdotDst{o}, 1)
        connect(m, 'Propellant_Flow', o, mdotDst{o}{k, 1}, mdotDst{o}{k, 2});
    end
end

% GROUND -----------------------------------------------------------------
gpos = get_param([m '/Sum4'], 'Position');
add_block('built-in/Subsystem', [m '/Ground_Contact'], 'Position', [gpos(1) - 260, gpos(2) + 80, gpos(1) - 120, gpos(2) + 200]);
gnd = {{'Inertial Position'; 'Inertial Velocity'; 'Inertial to Body'; 'Angular Velocity'}; {'Body Force'; 'Body Moment'}};
makeVariant(m, 'Ground_Contact', 'GROUND', {
    'Ground_None', @(c) zeroChoice(c, gnd{1}, gnd{2}, {'[0 0 0]', '[0 0 0]'})
    'Ground_Strut', @(c) refChoice(c, gnd, 'Ground_Contact')
    });
sixdofOut = [2 1 5 7];          % inertial position, velocity, DCM, body rates
for k = 1:4
    connect(m, '6DOF (Quaternion)', sixdofOut(k), 'Ground_Contact', k);
end
connect(m, 'Ground_Contact', 1, 'Sum4', 2);   % body force
connect(m, 'Ground_Contact', 2, 'Sum5', 2);   % body moment

set_param(m, 'InitFcn', '');   % FuelProps comes from buildSimInput when the fluid network is active
save_system(m);
close_system(m, 0);
fprintf('Built %s\n', out);
end

% ------------------------------------------------------------------------
function makeVariant(m, name, axis, choices)
% Turn block m/name (or an empty placeholder) into a Variant Subsystem whose
% choices are built by the given functions. Wiring of m/name is kept.
blk = [m '/' name];
A = variantAxes();
if size(choices, 1) ~= numel(A.(axis))
    error('build_hopper_6dof:choices', '%s has %d choices in variantAxes, %d here.', axis, numel(A.(axis)), size(choices, 1));
end
[ins, outs] = recordLines(blk);
pos = get_param(blk, 'Position');
tmp = [m '/' name '__choices'];
add_block('built-in/Subsystem', tmp);
for k = 1:size(choices, 1)
    choices{k, 2}([tmp '/' choices{k, 1}]);
end
p = portNames([tmp '/' choices{1, 1}]);
delete_block(blk);
add_block('simulink/Ports & Subsystems/Variant Subsystem', blk, 'Position', pos);
Simulink.SubSystem.deleteContents(blk);
for k = 1:numel(p.in)
    add_block('built-in/Inport', [blk '/' p.in{k}], 'Position', [20, 40 * k, 50, 40 * k + 14]);
end
for k = 1:numel(p.out)
    add_block('built-in/Outport', [blk '/' p.out{k}], 'Position', [500, 40 * k, 530, 40 * k + 14]);
end
for k = 1:size(choices, 1)
    c = [blk '/' choices{k, 1}];
    add_block([tmp '/' choices{k, 1}], c, 'Position', [200, 60 + 120 * (k - 1), 360, 160 + 120 * (k - 1)]);
    set_param(c, 'VariantControl', sprintf('VARIANT_%s == %d', axis, k));
end
delete_block(tmp);
set_param(blk, 'VariantControlMode', 'expression', 'VariantActivationTime', 'update diagram', ...
    'PropagateVariantConditions', 'off', 'AllowZeroVariantControls', 'off');
redrawLines(blk, ins, outs);
end

function c = lutChoice(c, props)
emptyChoice(c, props{1}, props{2});
add_block('simulink/Ports & Subsystems/Subsystem Reference', [c '/Thrust_to_Mdot_LUT'], ...
    'ReferencedSubsystem', 'Thrust_to_Mdot_LUT', 'Position', [200 80 340 160]);
add_line(c, 'Thrust_Cmd/1', 'Thrust/1');
add_line(c, 'Thrust_Cmd/1', 'Thrust_to_Mdot_LUT/1');
add_line(c, 'Thrust_to_Mdot_LUT/1', 'Ox_mdot/1');   % the LUT's first output is the ox mdot in v2's wiring
add_line(c, 'Thrust_to_Mdot_LUT/2', 'Fu_mdot/1');
end

function c = fluidChoice(c, props, src)
emptyChoice(c, props{1}, props{2});
add_block(src, [c '/Fluid_Network'], 'Position', [200 60 360 200]);
add_line(c, 'Thrust_Cmd/1', 'Fluid_Network/1');
add_line(c, 'Fluid_Network/1', 'Thrust/1');
add_line(c, 'Fluid_Network/2', 'Ox_mdot/1');
add_line(c, 'Fluid_Network/3', 'Fu_mdot/1');
end

function c = refChoice(c, ports, file)
emptyChoice(c, ports{1}, ports{2});
add_block('simulink/Ports & Subsystems/Subsystem Reference', [c '/' file], ...
    'ReferencedSubsystem', file, 'Position', [200 40 360 40 + 40 * max(numel(ports{1}), numel(ports{2}))]);
for k = 1:numel(ports{1}), add_line(c, sprintf('%s/1', ports{1}{k}), sprintf('%s/%d', file, k)); end
for k = 1:numel(ports{2}), add_line(c, sprintf('%s/%d', file, k), sprintf('%s/1', ports{2}{k})); end
end

function c = passChoice(c, ports)
emptyChoice(c, ports{1}, ports{2});
for k = 1:numel(ports{1})
    add_line(c, [ports{1}{k} '/1'], [ports{2}{k} '/1']);
end
end

function c = zeroChoice(c, ins, outs, values)
% Inputs end in terminators; each output is a constant.
emptyChoice(c, ins, outs);
for k = 1:numel(ins)
    t = sprintf('%s/Unused%d', c, k);
    add_block('built-in/Terminator', t, 'Position', [120, 40 * k, 140, 40 * k + 14]);
    add_line(c, [ins{k} '/1'], sprintf('Unused%d/1', k));
end
for k = 1:numel(outs)
    z = sprintf('Zero%d', k);
    add_block('built-in/Constant', [c '/' z], 'Value', values{k}, 'Position', [300, 40 * k, 360, 40 * k + 14]);
    add_line(c, [z '/1'], [outs{k} '/1']);
end
end

function emptyChoice(c, ins, outs)
if getSimulinkBlockHandle(c) < 0, add_block('built-in/Subsystem', c); end
for k = 1:numel(ins)
    add_block('built-in/Inport', [c '/' ins{k}], 'Position', [20, 40 * k, 50, 40 * k + 14]);
end
for k = 1:numel(outs)
    add_block('built-in/Outport', [c '/' outs{k}], 'Position', [500, 40 * k, 530, 40 * k + 14]);
end
end

function p = portNames(blk)
% Port block names in port order, for a subsystem or a subsystem reference
ins = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Inport');
outs = find_system(blk, 'SearchDepth', 1, 'BlockType', 'Outport');
p.in = sortByPort(ins);
p.out = sortByPort(outs);
end

function names = sortByPort(blks)
n = cellfun(@(b) str2double(get_param(b, 'Port')), blks);
[~, i] = sort(n);
names = cellfun(@(b) get_param(b, 'Name'), blks(i), 'uni', 0)';
end

function [ins, outs] = recordLines(blk)
% ins{p} = source port handle of input p; outs{p} = {dst port handles, line name}
ph = get_param(blk, 'PortHandles');
ins = cell(1, numel(ph.Inport)); outs = cell(1, numel(ph.Outport));
sys = get_param(blk, 'Parent');
for p = 1:numel(ph.Inport)
    l = get_param(ph.Inport(p), 'Line');
    if l > 0 && get_param(l, 'SrcPortHandle') > 0
        ins{p} = get_param(l, 'SrcPortHandle');
        delete_line(sys, ins{p}, ph.Inport(p));
    end
end
for p = 1:numel(ph.Outport)
    l = get_param(ph.Outport(p), 'Line');
    if l > 0
        d = get_param(l, 'DstPortHandle');
        outs{p} = {d(d > 0), get_param(l, 'Name')};
        delete_line(l);
    end
end
end

function redrawLines(blk, ins, outs)
sys = get_param(blk, 'Parent');
ph = get_param(blk, 'PortHandles');
for p = 1:numel(ins)
    if ~isempty(ins{p}), add_line(sys, ins{p}, ph.Inport(p), 'autorouting', 'smart'); end
end
for p = 1:numel(outs)
    if isempty(outs{p}), continue; end
    for d = outs{p}{1}(:)'
        l = add_line(sys, ph.Outport(p), d, 'autorouting', 'smart');
        set_param(l, 'Name', outs{p}{2});
    end
end
end

function deleteWithLines(blk)
% Inputs: remove only the branch into blk, so a source that also feeds
% other blocks keeps those lines. Outputs: the whole line is blk's own.
recordLines(blk);
delete_block(blk);
end

function disconnectInput(blk, port)
lh = get_param(blk, 'LineHandles');
h = lh.Inport(port);
if h > 0
    ph = get_param(blk, 'PortHandles');
    delete_line(get_param(blk, 'Parent'), get_param(h, 'SrcPortHandle'), ph.Inport(port));
end
end

function dst = destinations(blk, port)
% {dstName, dstPort; ...} for outport port of blk
ph = get_param(blk, 'PortHandles');
l = get_param(ph.Outport(port), 'Line');
dst = cell(0, 2);
if l < 0, return; end
d = get_param(l, 'DstBlockHandle'); dp = get_param(l, 'DstPortHandle');
for i = find(d > 0)'
    dst(end+1, :) = {get_param(d(i), 'Name'), get_param(dp(i), 'PortNumber')}; %#ok<AGROW>
end
end

function connect(m, srcName, srcPort, dstName, dstPort)
sp = get_param([m '/' strrep(srcName, '/', '//')], 'PortHandles');
dp = get_param([m '/' strrep(dstName, '/', '//')], 'PortHandles');
add_line(m, sp.Outport(srcPort), dp.Inport(dstPort), 'autorouting', 'smart');
end

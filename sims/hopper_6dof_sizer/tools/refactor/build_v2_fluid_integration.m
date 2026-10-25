function build_v2_fluid_integration()
%BUILD_V2_FLUID_INTEGRATION  Rebuild PR #26's fluid-network model on top of v2.
%   PR #26 enabled the Fluid Network2 subsystem inside hopper_6dof_NED_v3, a
%   v3 (error-quaternion LQR) model that needs sim_setup_v3. This script
%   applies the same change to the current master, hopper_6dof_NED_v2, so the
%   result runs from sim_setup:
%     - replace v2's commented-out Fluid Network2 with #26's version (enabled)
%     - feed it the thrust command from Subsystem2, as #26 did
%     - Actual Thrust drives the actuators, IC1 and the thrust logging
%     - Actual Ox/Fuel mdot replace the thrust->mdot lookup (Subsystem),
%       which is commented out
%     - #27's ground contact (strut spring/damper per leg) feeds body force
%       and moment into Sum4/Sum5, so the vehicle sits on the pad while the
%       fluid network spools up (~1 s to build thrust). v2's own "Ground "
%       block only drives scopes and stays commented out.
%   Writes models/hopper_6dof_NED_v2_fluid_integration.slx. Run with the
%   HopperSim project open.

src = 'hopper_6dof_NED_v3_fluid_integration';   % archive/models, PR #26
gnd = 'hopper_6dof_NED_v2_ground';              % models/, PR #27
base = 'hopper_6dof_NED_v2';
new = 'hopper_6dof_NED_v2_fluid_integration';
out = hopperRoot('models', [new '.slx']);

addpath(hopperRoot('archive', 'models'));
cleanup = onCleanup(@() rmpath(hopperRoot('archive', 'models')));
load_system(src);
load_system(base);
if bdIsLoaded(new), close_system(new, 0); end
save_system(base, out);          % v2 copy under the new name
close_system(base, 0);
m = new;

% Swap in #26's fluid network at the same position
fn = [m '/Fluid Network2'];
pos = get_param(fn, 'Position');
deleteLines(fn);
delete_block(fn);
add_block([src '/Fluid Network2'], fn, 'Position', pos);
set_param(fn, 'Commented', 'off');

% The thrust -> mdot lookup is replaced by the fluid network
lut = [m '/Subsystem'];
lutDst = destinations(lut);      % where its Ox/Fuel mdot outputs went
deleteLines(lut);
set_param(lut, 'Commented', 'on');

% Thrust command now goes to the fluid network instead of the actuators
thrustDst = {'IC1', 1; 'Control Inputs -> Fxyz & Mxyz', 1; 'To Workspace7', 1; 'Thrust', 1};
for k = 1:size(thrustDst, 1)
    disconnectInput([m '/' thrustDst{k, 1}], thrustDst{k, 2});
end
connect(m, 'Subsystem2', 1, 'Fluid Network2', 1);

% Actual Thrust (out 1) and Actual Ox/Fuel mdot (outs 2/3)
for k = 1:size(thrustDst, 1)
    connect(m, 'Fluid Network2', 1, thrustDst{k, 1}, thrustDst{k, 2});
end
for o = 1:2
    for d = lutDst{o}'
        connect(m, 'Fluid Network2', o + 1, d{1}, d{2});
    end
end

% Ground contact from #27: same 6DOF inputs, forces/moments into Sum4/Sum5
load_system(gnd);
gpos = get_param([gnd '/Ground'], 'Position');
add_block([gnd '/Ground'], [m '/Ground_Contact'], 'Position', gpos + [0 150 0 150]);
close_system(gnd, 0);
sixdofOut = [2 1 5 7];          % Inertial position, velocity, DCM, body rates
for k = 1:4
    connect(m, '6DOF (Quaternion)', sixdofOut(k), 'Ground_Contact', k);
end
connect(m, 'Ground_Contact', 1, 'Sum4', 2);   % body force
connect(m, 'Ground_Contact', 2, 'Sum5', 2);   % body moment

% The fluid network's tank blocks read FuelProps from FuelPropsData.mat (data/)
set_param(m, 'InitFcn', 'load(''FuelPropsData.mat'');');
save_system(m);
close_system(m, 0);
close_system(src, 0);
fprintf('Built %s\n', out);
end

function deleteLines(blk)
lh = get_param(blk, 'LineHandles');
for h = [lh.Inport, lh.Outport]
    if h > 0, delete_line(h); end
end
end

function disconnectInput(blk, port)
lh = get_param(blk, 'LineHandles');
h = lh.Inport(port);
if h > 0
    % a branched line: remove only the branch into this port
    ph = get_param(blk, 'PortHandles');
    delete_line(get_param(blk, 'Parent'), get_param(h, 'SrcPortHandle'), ph.Inport(port));
end
end

function dst = destinations(blk)
% dst{o} = {dstName, dstPort; ...} for each outport o of blk
pc = get_param(blk, 'PortConnectivity');
pc = pc(arrayfun(@(p) isempty(p.SrcBlock), pc));
dst = cell(numel(pc), 1);
for o = 1:numel(pc)
    dst{o} = [arrayfun(@(h) get_param(h, 'Name'), pc(o).DstBlock(:), 'uni', 0), num2cell(pc(o).DstPort(:) + 1)];
end
end

function connect(m, srcName, srcPort, dstName, dstPort)
% Port handles, because block names like 'CG / MOI Script v3' contain '/'
sp = get_param([m '/' strrep(srcName, '/', '//')], 'PortHandles');
dp = get_param([m '/' strrep(dstName, '/', '//')], 'PortHandles');
add_line(m, sp.Outport(srcPort), dp.Inport(dstPort), 'autorouting', 'smart');
end

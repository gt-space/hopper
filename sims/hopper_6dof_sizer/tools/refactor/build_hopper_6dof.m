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
%     Controller       CONTROLLER  lqr_sched: LQR_Gain_Scheduler with v2's
%                                       clock and gain tables
%                                  lqi: LQI-DEV's regulator from
%                                       models/dev/hopper_6dof_NED_v4.slx
%                                  lqr_sched_legacy: the March scheduler code
%     Sensors          FDIR on or  the shared sensor model + sensor faults
%                      NAVIGATION  (src/sensors/sensor_frame.m) at the IMU
%                                  rate; otherwise zeros
%     Navigation       NAVIGATION  none: truth to the controller and FDIR
%                                  gmekf_shadow: the nav team's Hopper_GMEKF
%                                       (Model block) on the Sensors block,
%                                       logged; truth still flies
%                                  gmekf: the controller and FDIR fly on
%                                       the estimate
%     FDIR             FDIR        off / on (src/fdir/fdir_frame.m)
%     Actuator_Faults  FAULT_INJECTION  off / on (src/faults/actuator_gate.m)
%
%   VEHICLE has no block: it picks the parameter set (hopperWorkspace).
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

% CONTROLLER: between the 13-state vector (Mux) and the actuator split --
sched = [m '/LQR_Gain_Scheduler'];
[schedIn, schedOut] = recordLines(sched);
gainSrc = cellfun(@(h) get_param(get_param(h, 'Parent'), 'Name'), schedIn([1 3:7]), 'uni', 0);  % Clock, tgrid, xref, K1flat, K2grid, unom
xSrc = schedIn{2};
cpos = get_param(sched, 'Position');
delete_block(sched);
add_block('built-in/Subsystem', [m '/Controller'], 'Position', cpos);
v4 = 'hopper_6dof_NED_v4';
load_system(v4);
closeV4 = onCleanup(@() close_system(v4, 0));
makeVariant(m, 'Controller', 'CONTROLLER', {
    'Controller_LQR_Sched', @(c) schedChoice(c, m, gainSrc, 'current')
    'Controller_LQI', @(c) lqiChoice(c, v4)
    'Controller_LQR_Sched_Legacy', @(c) schedChoice(c, m, gainSrc, 'legacy')
    });
for g = gainSrc   % the gain tables and clock now live in the scheduler choices
    lh = get_param([m '/' g{1}], 'LineHandles');
    if all(lh.Outport < 0), delete_block([m '/' g{1}]); end   % keep one that feeds something else
end
ph = get_param([m '/Controller'], 'PortHandles');

% SENSORS and NAVIGATION: truth -> shared Sensors block (500 Hz) -> nav
% filter -> the state the controller and FDIR fly on
add_block('built-in/Subsystem', [m '/Sensors'], 'Position', cpos + [-260 140 -160 220]);
sensPorts = {{'x'; 'f_b'; 'dcm_be'; 'cg'}; {'meas'; 'valid'; 'stamp'}};
sensOn = 'VARIANT_FDIR == 2 || VARIANT_NAVIGATION > 1';
makeVariant(m, 'Sensors', '', {
    'Sensors_Off', @(c) zeroChoice(c, sensPorts{1}, sensPorts{2}, {'zeros(26,1)', 'zeros(9,1)', 'zeros(9,1)'})
    'Sensors_On',  @(c) sensorsChoice(c, sensPorts)
    }, {['~(' sensOn ')'], sensOn});
sensSrc = {xSrc, port(m, 'Slosh Specific Force', 1), port(m, '6DOF (Quaternion)', 5), port(m, 'Mass_Properties', 1)};
sph = get_param([m '/Sensors'], 'PortHandles');
for k = 1:numel(sensSrc)
    add_line(m, sensSrc{k}, sph.Inport(k), 'autorouting', 'smart');
end
add_block('built-in/Subsystem', [m '/Navigation'], 'Position', cpos + [-130 0 -30 120]);
navPorts = {{'x_true'; 'dcm_true'; 'meas'; 'valid'; 'stamp'; 'cg'; 'mass'; 'MOI'; 'Mxyz'}; {'x'; 'dcm'; 'xn'; 'P'}};
makeVariant(m, 'Navigation', 'NAVIGATION', {
    'Navigation_None',   @(c) passChoiceFirst2(c, navPorts)
    'Navigation_Shadow', @(c) navChoice(c, navPorts, false)
    'Navigation_GMEKF',  @(c) navChoice(c, navPorts, true)
    });
navSrc = {xSrc, port(m, '6DOF (Quaternion)', 5), sph.Outport(1), sph.Outport(2), sph.Outport(3), ...
    port(m, 'Mass_Properties', 1), port(m, 'Saturation1', 1), port(m, 'Mass_Properties', 2), port(m, 'Actuators_TVC_RCS', 2)};
nph = get_param([m '/Navigation'], 'PortHandles');
for k = 1:numel(navSrc)
    add_line(m, navSrc{k}, nph.Inport(k), 'autorouting', 'smart');
end
add_line(m, nph.Outport(1), ph.Inport(1), 'autorouting', 'smart');
for d = schedOut{1}{1}(:)'
    add_line(m, ph.Outport(1), d, 'autorouting', 'smart');
end

% FDIR and FAULT_INJECTION: controller -> FDIR -> actuator gates -> actuators
uLine = get_param(ph.Outport(1), 'Line');
uDst = get_param(uLine, 'DstPortHandle');
delete_line(uLine);
cpos = get_param([m '/Controller'], 'Position');
add_block('built-in/Subsystem', [m '/FDIR'], 'Position', cpos + [240 0 300 80]);
fdirPorts = {{'u_nom'; 'x'; 'f_b'; 'dcm_be'; 'thrust'; 'tvc_pitch'; 'tvc_yaw'; 'Mxyz'; 'mass'; 'meas'; 'valid'; 'stamp'; 'xn'; 'P'; 'cg'}; {'u'}};
makeVariant(m, 'FDIR', 'FDIR', {
    'FDIR_Off', @(c) passFirst(c, fdirPorts)
    'FDIR_On',  @(c) fdirChoice(c, fdirPorts)
    });
add_block('built-in/Subsystem', [m '/Actuator_Faults'], 'Position', cpos + [440 0 520 40]);
gatePorts = {{'u_in'}; {'u'}};
makeVariant(m, 'Actuator_Faults', 'FAULT_INJECTION', {
    'Actuator_Faults_Off', @(c) passFirst(c, gatePorts)
    'Actuator_Faults_On',  @(c) actuatorGateChoice(c, gatePorts)
    });
fdirSrc = {ph.Outport(1), nph.Outport(1), port(m, 'Slosh Specific Force', 1), nph.Outport(2), ...
    port(m, 'Propulsion', 1), port(m, 'Actuators_TVC_RCS', 3), port(m, 'Actuators_TVC_RCS', 4), ...
    port(m, 'Actuators_TVC_RCS', 2), port(m, 'Saturation1', 1), sph.Outport(1), sph.Outport(2), sph.Outport(3), ...
    nph.Outport(3), nph.Outport(4), port(m, 'Mass_Properties', 1)};
fph = get_param([m '/FDIR'], 'PortHandles');
for k = 1:numel(fdirSrc)
    add_line(m, fdirSrc{k}, fph.Inport(k), 'autorouting', 'smart');
end
gph = get_param([m '/Actuator_Faults'], 'PortHandles');
add_line(m, fph.Outport(1), gph.Inport(1), 'autorouting', 'smart');
for d = uDst(uDst > 0)'
    add_line(m, gph.Outport(1), d, 'autorouting', 'smart');
end

linkNavigationDictionary(m);   % the nav filter's bus types
set_param(m, 'InitFcn', '');   % FuelProps comes from buildSimInput when the fluid network is active
save_system(m);
close_system(m, 0);
fprintf('Built %s\n', out);
end

% ------------------------------------------------------------------------
function makeVariant(m, name, axis, choices, conds)
% Turn block m/name (or an empty placeholder) into a Variant Subsystem whose
% choices are built by the given functions. Wiring of m/name is kept.
% Choice k is active when VARIANT_<axis> == k, or when conds{k} holds if
% conditions are given (a block that depends on several axes).
blk = [m '/' name];
A = variantAxes();
if nargin < 5
    if size(choices, 1) ~= numel(A.(axis))
        error('build_hopper_6dof:choices', '%s has %d choices in variantAxes, %d here.', axis, numel(A.(axis)), size(choices, 1));
    end
    conds = arrayfun(@(k) sprintf('VARIANT_%s == %d', axis, k), 1:size(choices, 1), 'uni', 0);
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
    set_param(c, 'VariantControl', conds{k});
end
delete_block(tmp);
set_param(blk, 'VariantControlMode', 'expression', 'VariantActivationTime', 'update diagram', ...
    'PropagateVariantConditions', 'off', 'AllowZeroVariantControls', 'off');
redrawLines(blk, ins, outs);
end

function c = schedChoice(c, m, gainSrc, kind)
% Time-scheduled LQR: clock and gain tables from v2's top level, and the
% scheduler (the shared file, or the March code inline for 'legacy')
emptyChoice(c, {'x'}, {'u'});
for g = gainSrc
    add_block([m '/' g{1}], [c '/' g{1}]);
end
switch kind
    case 'current'
        s = 'LQR_Gain_Scheduler';
        add_block('simulink/Ports & Subsystems/Subsystem Reference', [c '/' s], ...
            'ReferencedSubsystem', s, 'Position', [250 40 400 240]);
    case 'legacy'
        s = 'Gain_Selector_Legacy';
        add_block('simulink/User-Defined Functions/MATLAB Function', [c '/' s], 'Position', [250 40 400 240]);
        ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', [c '/' s]);
        ch.Script = legacyGainSelector();
end
% scheduler inputs: t, x, tgrid, xref, K1flat, K2grid, unom
add_line(c, [gainSrc{1} '/1'], [s '/1']);
add_line(c, 'x/1', [s '/2']);
for k = 2:numel(gainSrc)
    add_line(c, [gainSrc{k} '/1'], sprintf('%s/%d', s, k + 1));
end
add_line(c, [s '/1'], 'u/1');
end

function c = lqiChoice(c, v4)
% LQI-DEV's regulator (models/dev/hopper_6dof_NED_v4.slx):
%   u = -Klqi * [pos - cmd; vel; rates; MRP; integral(pos - cmd)]
% with the attitude as MRPs (v4's 'q to MRP' and 'MRP Correction') and a
% zero XYZ command. v4's own 'Apply Gains' block computes these terms but
% joins them with Mux blocks (16 elements into the 4 actuator inputs), so
% it doesn't compile; the law is built here directly with one gain block.
emptyChoice(c, {'x'}, {'u'});
add_block('built-in/Selector', [c '/Position'], 'InputPortWidth', '13', 'Indices', '1:3', 'Position', [100 30 140 60]);
add_block('built-in/Selector', [c '/Vel_Rates'], 'InputPortWidth', '13', 'Indices', '4:9', 'Position', [100 90 140 120]);
add_block('built-in/Selector', [c '/Quaternion'], 'InputPortWidth', '13', 'Indices', '10:13', 'Position', [100 150 140 180]);
for b = {'q to MRP', 'MRP Correction'}
    add_block([v4 '/' b{1}], [c '/' b{1}]);
end
add_block('built-in/Constant', [c '/Position_Command'], 'Value', '[0 0 0]', 'Position', [100 230 160 250]);
add_block('built-in/Sum', [c '/Position_Error'], 'Inputs', '+-', 'Position', [220 30 250 60]);
add_block('built-in/Integrator', [c '/Position_Error_Integral'], 'InitialCondition', '[0 0 0]', 'Position', [300 230 330 260]);
add_block('built-in/Mux', [c '/LQI_State'], 'Inputs', '4', 'Position', [400 30 405 260]);
add_block('built-in/Gain', [c '/LQI_Gain'], 'Gain', '-Klqi', 'Multiplication', 'Matrix(K*u)', 'Position', [450 125 510 165]);
add_line(c, 'x/1', 'Position/1');
add_line(c, 'x/1', 'Vel_Rates/1');
add_line(c, 'x/1', 'Quaternion/1');
add_line(c, 'Quaternion/1', 'q to MRP/1');
add_line(c, 'q to MRP/1', 'MRP Correction/1');
add_line(c, 'Position/1', 'Position_Error/1');
add_line(c, 'Position_Command/1', 'Position_Error/2');
add_line(c, 'Position_Error/1', 'Position_Error_Integral/1');
add_line(c, 'Position_Error/1', 'LQI_State/1');
add_line(c, 'Vel_Rates/1', 'LQI_State/2');
add_line(c, 'MRP Correction/1', 'LQI_State/3');
add_line(c, 'Position_Error_Integral/1', 'LQI_State/4');
add_line(c, 'LQI_State/1', 'LQI_Gain/1');
add_line(c, 'LQI_Gain/1', 'u/1');
end

function c = passFirst(c, ports)
% first input straight to the output; the others end in terminators
emptyChoice(c, ports{1}, ports{2});
add_line(c, [ports{1}{1} '/1'], [ports{2}{1} '/1']);
for k = 2:numel(ports{1})
    add_block('built-in/Terminator', sprintf('%s/Unused%d', c, k), 'Position', [120, 40 * k, 140, 40 * k + 14]);
    add_line(c, [ports{1}{k} '/1'], sprintf('Unused%d/1', k));
end
end

function c = fdirChoice(c, ports)
% Monitors, executive and mode manager at 100 Hz (fdir_frame) on the
% shared Sensors block's measurements (rate transitions from its 500 Hz),
% the R1 low-bandwidth scheduler, and the command select.
% The frame's outputs reach the command one frame later (Unit Delays): a
% mode change takes effect next frame (GNC PDR slide 78), and it keeps the
% actuator feedback FDIR reads out of an algebraic loop.
emptyChoice(c, ports{1}, ports{2});
fr = 'fdir_params.S.frame';
add_block('built-in/Clock', [c '/Clock'], 'Position', [100 520 120 540]);
rt = {'meas', 'valid', 'stamp', 'xn', 'P'};
for k = 1:numel(rt)
    add_block('built-in/RateTransition', [c '/RT_' rt{k}], 'OutPortSampleTime', fr, ...
        'Position', [180, 600 + 40 * k, 210, 620 + 40 * k]);
    add_line(c, [rt{k} '/1'], ['RT_' rt{k} '/1']);
end

f = [c '/FDIR_Frame'];
add_block('simulink/User-Defined Functions/MATLAB Function', f, 'Position', [300 40 460 460]);
setScript(f, {
    'function [abort, safed, degrade, x_ref, mode, level, first_fault, status, raw, nav_abort, nav_nis] = fcn(t, x, f_b, dcm_be, thrust, tvc_pitch, tvc_yaw, Mxyz, mass, u_nom, u_cmd, meas, valid, stamp, xn, P, cg, fdir_params)'
    '% One 100 Hz FDIR frame (src/fdir/fdir_frame.m, docs/FDIR_ARCHITECTURE.md)'
    'persistent st'
    'if isempty(st)'
    '    st = fdir_frame_state(size(fdir_params.cat, 1), size(fdir_params.faults, 1), fdir_params.S.frame);'
    'end'
    'in.t = t; in.x = x(:); in.f_b = f_b(:); in.dcm_be = dcm_be;'
    'in.meas = meas(:); in.valid = valid(:); in.stamp = stamp(:);'
    'in.xn = xn(:); in.P = P; in.cg = cg(:);'
    'in.thrust = thrust; in.gimbal_fb = [tvc_pitch; tvc_yaw]; in.rcs_torque = Mxyz(1);'
    'in.mass = mass; in.u_nom = u_nom(:); in.u_cmd = u_cmd(:);'
    '[y, st] = fdir_frame(in, fdir_params, st);'
    'abort = double(y.abort); safed = double(y.safed); degrade = double(y.degrade);'
    'x_ref = y.x_ref; mode = y.mode; level = y.level; first_fault = double(y.first_fault);'
    'status = y.status; raw = y.raw; nav_abort = double(y.nav_abort); nav_nis = y.nav_nis;'
    'end'}, {'fdir_params'}, '0.01');   % the 100 Hz FDIR frame (sim_sensor_params.frame)
% FDIR_Frame inputs: t, then x f_b dcm_be thrust tvc_pitch tvc_yaw Mxyz mass
% (ports 2-9), u_nom (10), u_cmd (11, below), meas valid stamp (12-14)
add_line(c, 'Clock/1', 'FDIR_Frame/1');
for k = 2:9
    add_line(c, [ports{1}{k} '/1'], sprintf('FDIR_Frame/%d', k));
end
add_line(c, 'u_nom/1', 'FDIR_Frame/10');
add_line(c, 'RT_meas/1', 'FDIR_Frame/12');
add_line(c, 'RT_valid/1', 'FDIR_Frame/13');
add_line(c, 'RT_stamp/1', 'FDIR_Frame/14');
add_line(c, 'RT_xn/1', 'FDIR_Frame/15');
add_line(c, 'RT_P/1', 'FDIR_Frame/16');
add_line(c, 'cg/1', 'FDIR_Frame/17');

% R1: the same scheduler with the low-bandwidth gains
lo = [c '/LQR_Low_Bandwidth'];
add_block('simulink/Ports & Subsystems/Subsystem Reference', lo, 'ReferencedSubsystem', 'LQR_Gain_Scheduler', ...
    'Position', [300 520 450 700]);
add_line(c, 'Clock/1', 'LQR_Low_Bandwidth/1');
add_line(c, 'x/1', 'LQR_Low_Bandwidth/2');
tables = {'tgrid', 'Target_Trajectory', 'K1flat_lo', 'K2grid', 'unom'};
for k = 1:numel(tables)
    b = ['Table_' tables{k}];
    add_block('built-in/Constant', [c '/' b], 'Value', tables{k}, 'Position', [180, 500 + 30 * k, 240, 515 + 30 * k]);
    add_line(c, [b '/1'], sprintf('LQR_Low_Bandwidth/%d', k + 2));
end

% one-frame delays into the command select
outs = {'abort', 'safed', 'degrade', 'x_ref', 'nav_abort'};
for k = 1:numel(outs)
    d = ['Delay_' outs{k}];
    add_block('built-in/UnitDelay', [c '/' d], 'SampleTime', fr, 'InitialCondition', '0', ...
        'Position', [520, 40 + 50 * k, 550, 70 + 50 * k]);
    add_line(c, sprintf('FDIR_Frame/%d', find(strcmp({'abort', 'safed', 'degrade', 'x_ref', '', '', '', '', '', 'nav_abort'}, outs{k}))), [d '/1']);
end
set_param([c '/Delay_x_ref'], 'InitialCondition', 'zeros(13,1)');

s = [c '/Command_Select'];
add_block('simulink/User-Defined Functions/MATLAB Function', s, 'Position', [640 40 780 320]);
setScript(s, {
    'function u = fcn(u_nom, u_lo, x, x_ref, mass, abort, safed, degrade, nav_abort, K_abort)'
    '% The response ladder applied to the command (src/fdir/fdir_command_select.m)'
    'mode.abort = abort > 0.5; mode.safed = safed > 0.5; mode.degrade = degrade > 0.5; mode.nav_abort = nav_abort > 0.5;'
    'u = fdir_command_select(u_nom, u_lo, x, x_ref, mass, mode, K_abort, 9.80665);'
    'end'}, {'K_abort'}, '-1');
add_line(c, 'u_nom/1', 'Command_Select/1');
add_line(c, 'LQR_Low_Bandwidth/1', 'Command_Select/2');
add_line(c, 'x/1', 'Command_Select/3');
add_line(c, 'Delay_x_ref/1', 'Command_Select/4');
add_line(c, 'mass/1', 'Command_Select/5');
add_line(c, 'Delay_abort/1', 'Command_Select/6');
add_line(c, 'Delay_safed/1', 'Command_Select/7');
add_line(c, 'Delay_degrade/1', 'Command_Select/8');
add_line(c, 'Delay_nav_abort/1', 'Command_Select/9');
add_line(c, 'Command_Select/1', 'u/1');
% the command actually sent, back to the monitors one frame later
add_block('built-in/UnitDelay', [c '/Delay_u_cmd'], 'SampleTime', fr, 'InitialCondition', 'zeros(4,1)', 'Position', [820 360 850 390]);
add_line(c, 'Command_Select/1', 'Delay_u_cmd/1');
add_line(c, 'Delay_u_cmd/1', 'FDIR_Frame/11');

% logged for runs and Monte Carlo summaries
logs = {1, 'fdir_abort'; 2, 'fdir_safed'; 3, 'fdir_degrade'; 5, 'fdir_mode'; 6, 'fdir_level'; ...
        7, 'fdir_first_fault'; 8, 'fdir_status'; 9, 'fdir_raw'; 11, 'fdir_nav_nis'};
for k = 1:size(logs, 1)
    b = ['Log_' logs{k, 2}];
    add_block('built-in/ToWorkspace', [c '/' b], 'VariableName', logs{k, 2}, 'SaveFormat', 'Timeseries', 'MaxDataPoints', 'inf', ...
        'Position', [880, 40 + 40 * k, 960, 60 + 40 * k]);
    add_line(c, sprintf('FDIR_Frame/%d', logs{k, 1}), [b '/1']);
end
end

function c = actuatorGateChoice(c, ports)
% Actuator faults from the fault table. The gate that touches the command
% is stateless and continuous, so with no fault active the flight is
% untouched; its memory (command at onset, slow-mode rate limit) is a
% separate block sampled at the solver step. A sampled command changed the
% thrust path's If/Merge cutoff logic even with no fault.
emptyChoice(c, ports{1}, ports{2});
add_block('built-in/Clock', [c '/Clock'], 'Position', [100 160 120 180]);
s = [c '/Actuator_Gate_State'];
add_block('simulink/User-Defined Functions/MATLAB Function', s, 'Position', [220 120 360 220]);
setScript(s, {
    'function [hold, slow] = fcn(u_in, t, fdir_params)'
    '% Memory for the actuator gate (src/faults/actuator_gate_state.m)'
    'persistent g'
    'if isempty(g)'
    '    g = gate_state(size(fdir_params.faults, 1), 1e-3);'
    'end'
    '[hold, slow, g] = actuator_gate_state(u_in(:), t, fdir_params.faults, g);'
    'end'}, {'fdir_params'}, '1e-3');
f = [c '/Actuator_Gate'];
add_block('simulink/User-Defined Functions/MATLAB Function', f, 'Position', [420 40 560 200]);
setScript(f, {
    'function u = fcn(u_in, t, hold, slow, fdir_params)'
    '% Gimbal, throttle and RCS faults (src/faults/actuator_gate.m)'
    'u = u_in;'
    'if fdir_params.inject'
    '    u(:) = actuator_gate(u_in(:), t, fdir_params.faults, fdir_params.mon.gimbal_max, hold, slow);'
    'end'
    'end'}, {'fdir_params'}, '-1');
add_line(c, 'u_in/1', 'Actuator_Gate_State/1');
add_line(c, 'Clock/1', 'Actuator_Gate_State/2');
add_line(c, 'u_in/1', 'Actuator_Gate/1');
add_line(c, 'Clock/1', 'Actuator_Gate/2');
add_line(c, 'Actuator_Gate_State/1', 'Actuator_Gate/3');
add_line(c, 'Actuator_Gate_State/2', 'Actuator_Gate/4');
add_line(c, 'Actuator_Gate/1', 'u/1');
end

function c = sensorsChoice(c, ports)
% The shared Sensors block (src/sensors/sensor_frame.m): truth -> sensor
% model -> sensor faults, at the IMU rate. FDIR and the nav filter read it.
emptyChoice(c, ports{1}, ports{2});
add_block('built-in/Clock', [c '/Clock'], 'Position', [100 240 120 260]);
add_block('simulink/Sources/Random Number', [c '/Sensor_Noise'], 'Mean', 'zeros(26,1)', ...
    'Variance', 'ones(26,1)', 'Seed', 'sensor_params.S.seed + (0:25)''', 'SampleTime', 'sensor_params.dt', ...
    'Position', [100 300 140 330]);
f = [c '/Sensor_Frame'];
add_block('simulink/User-Defined Functions/MATLAB Function', f, 'Position', [240 40 400 330]);
setScript(f, {
    'function [meas, valid, stamp] = fcn(t, noise, x, f_b, dcm_be, cg, sensor_params)'
    '% Sensors and sensor faults, one step (src/sensors/sensor_frame.m)'
    'persistent st'
    'if isempty(st)'
    '    st = sensor_frame_state(size(sensor_params.faults, 1), sensor_params.dt);'
    'end'
    '[meas, v, stamp, ~, st] = sensor_frame(x(:), f_b(:), dcm_be, cg(:), noise(:), t, sensor_params, st);'
    'valid = double(v);'
    'end'}, {'sensor_params'}, 'sensor_params.dt');
add_line(c, 'Clock/1', 'Sensor_Frame/1');
add_line(c, 'Sensor_Noise/1', 'Sensor_Frame/2');
for k = 1:numel(ports{1})
    add_line(c, [ports{1}{k} '/1'], sprintf('Sensor_Frame/%d', k + 2));
end
for k = 1:numel(ports{2})
    add_line(c, sprintf('Sensor_Frame/%d', k), [ports{2}{k} '/1']);
end
end

function c = passChoiceFirst2(c, ports)
% Navigation none: truth state and DCM straight through; no estimate
emptyChoice(c, ports{1}, ports{2});
add_line(c, [ports{1}{1} '/1'], [ports{2}{1} '/1']);
add_line(c, [ports{1}{2} '/1'], [ports{2}{2} '/1']);
add_block('built-in/Constant', [c '/No_Estimate'], 'Value', 'zeros(20,1)', 'Position', [300 120 360 140]);
add_block('built-in/Constant', [c '/No_Covariance'], 'Value', 'zeros(19)', 'Position', [300 160 360 180]);
add_line(c, 'No_Estimate/1', [ports{2}{3} '/1']);
add_line(c, 'No_Covariance/1', [ports{2}{4} '/1']);
for k = 3:numel(ports{1})
    add_block('built-in/Terminator', sprintf('%s/Unused%d', c, k), 'Position', [120, 40 * k, 140, 40 * k + 14]);
    add_line(c, [ports{1}{k} '/1'], sprintf('Unused%d/1', k));
end
end

function c = navChoice(c, ports, closed)
% The nav team's GMEKF (subsystems/navigation/Filters/Hopper_GMEKF.slx, a
% Model block: their file, never copied) between our adapter blocks
% (src/navigation). closed: the controller and FDIR get the estimate;
% otherwise (shadow) truth, with the estimate logged beside it.
emptyChoice(c, ports{1}, ports{2});
add_block('built-in/UnitDelay', [c '/Delay_Mxyz'], 'SampleTime', 'nav_dt', 'InitialCondition', 'zeros(3,1)', ...
    'Position', [100 420 130 450]);   % the torque the flight computer commanded last step
add_line(c, 'Mxyz/1', 'Delay_Mxyz/1');
f = [c '/Nav_Inputs'];
add_block('simulink/User-Defined Functions/MATLAB Function', f, 'Position', [220 40 380 460]);
setScript(f, {
    'function [accel, angvel, baro, gps, lidar, dgps, lidarFrames, gpsFrame, massProps, torque] = fcn(meas, valid, stamp, cg, mass, MOI, Mxyz, nav_params)'
    '% The Sensors block as the filter''s buses (src/navigation/nav_measurements.m)'
    '[accel, angvel, baro, gps, lidar, dgps, ~, lidarFrames, gpsFrame, massProps, torque] = nav_measurements( ...'
    '    meas(:), valid(:) > 0.5, stamp(:), cg(:), mass, reshape(MOI, 3, 3), Mxyz(:), nav_params);'
    'end'}, {'nav_params'}, 'nav_dt');
busOut = {'accel', 'accelBus'; 'angvel', 'angvelBus'; 'baro', 'baroBus'; 'gps', 'gpsBus'; 'lidar', 'lidarBus'; ...
          'dgps', 'diffGpsBus'; 'lidarFrames', 'lidarFrameBus'; 'gpsFrame', 'gpsFrameBus'; 'massProps', 'massProperties'};
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', f);
for k = 1:size(busOut, 1)
    d = ch.find('-isa', 'Stateflow.Data', 'Name', busOut{k, 1});
    d.DataType = ['Bus: ' busOut{k, 2}];
end
for k = 3:7   % meas valid stamp cg mass
    add_line(c, [ports{1}{k} '/1'], sprintf('Nav_Inputs/%d', k - 2));
end
add_line(c, 'MOI/1', 'Nav_Inputs/6');
add_line(c, 'Delay_Mxyz/1', 'Nav_Inputs/7');

inherit = 'Inherit: Inherit from ''Constant value''';
consts = {'Init', 'nav_init', 'Bus: InitConditions'; 'Reset', 'nav_reset', inherit; ...
          'Step', 'nav_dt', inherit; 'Process_Noise', 'nav_Qc', inherit; 'Gravity', 'nav_g', inherit};
for k = 1:size(consts, 1)
    add_block('built-in/Constant', [c '/' consts{k, 1}], 'Value', consts{k, 2}, 'OutDataTypeStr', consts{k, 3}, ...
        'Position', [420, 500 + 30 * k, 480, 515 + 30 * k]);
end
add_block('simulink/Ports & Subsystems/Model', [c '/GMEKF'], 'ModelName', 'Hopper_GMEKF', 'Position', [560 40 720 560]);
% Hopper_GMEKF inputs: init reset dt Qc mass_properties gravity torque
% meas_posvel meas_accel meas_angvel meas_baro meas_LID lidar_frames
% meas_diffGPS gps_frame
wires = {'Init/1', 1; 'Reset/1', 2; 'Step/1', 3; 'Process_Noise/1', 4; 'Nav_Inputs/9', 5; 'Gravity/1', 6; ...
         'Nav_Inputs/10', 7; 'Nav_Inputs/4', 8; 'Nav_Inputs/1', 9; 'Nav_Inputs/2', 10; 'Nav_Inputs/3', 11; ...
         'Nav_Inputs/5', 12; 'Nav_Inputs/7', 13; 'Nav_Inputs/6', 14; 'Nav_Inputs/8', 15};
for k = 1:size(wires, 1)
    add_line(c, wires{k, 1}, sprintf('GMEKF/%d', wires{k, 2}));
end

o = [c '/Nav_Outputs'];
add_block('simulink/User-Defined Functions/MATLAB Function', o, 'Position', [800 40 940 200]);
setScript(o, {
    'function [x, dcm, xn, P_diag] = fcn(state_estimate, covariance, cg, nav_params)'
    '% The estimate in the controller''s convention (src/navigation/nav_state_to_sim.m)'
    'xn = state_estimate(:);'
    '[x, dcm] = nav_state_to_sim(xn, cg(:), nav_params.F);'
    'P_diag = diag(covariance);'
    'end'}, {'nav_params'}, '-1');
add_line(c, 'GMEKF/1', 'Nav_Outputs/1');
add_line(c, 'GMEKF/2', 'Nav_Outputs/2');
add_line(c, 'cg/1', 'Nav_Outputs/3');
add_line(c, 'Nav_Outputs/3', 'xn/1');     % the raw estimate and covariance, for FDIR's estimator checks
add_line(c, 'GMEKF/2', 'P/1');

add_block('built-in/Sum', [c '/Estimate_Error'], 'Inputs', '+-', 'Position', [1000 230 1030 260]);
add_line(c, 'Nav_Outputs/1', 'Estimate_Error/1');
add_line(c, 'x_true/1', 'Estimate_Error/2');
logs = {'Nav_Outputs/1', 'nav_estimate'; 'Estimate_Error/1', 'nav_error'; 'Nav_Outputs/3', 'nav_state'; 'Nav_Outputs/4', 'nav_P_diag'};
for k = 1:size(logs, 1)
    b = ['Log_' logs{k, 2}];
    add_block('built-in/ToWorkspace', [c '/' b], 'VariableName', logs{k, 2}, 'SaveFormat', 'Timeseries', 'MaxDataPoints', 'inf', ...
        'Position', [1080, 200 + 40 * k, 1160, 220 + 40 * k]);
    add_line(c, logs{k, 1}, [b '/1']);
end
if closed
    add_line(c, 'Nav_Outputs/1', 'x/1');
    add_line(c, 'Nav_Outputs/2', 'dcm/1');
    add_block('built-in/Terminator', [c '/Unused_dcm_true'], 'Position', [120 80 140 94]);
    add_line(c, 'dcm_true/1', 'Unused_dcm_true/1');
else
    add_line(c, 'x_true/1', 'x/1');
    add_line(c, 'dcm_true/1', 'dcm/1');
    add_block('built-in/Terminator', [c '/Unused_dcm'], 'Position', [980 100 1000 114]);
    add_line(c, 'Nav_Outputs/2', 'Unused_dcm/1');
end
end

function linkNavigationDictionary(m)
% The filter's buses live in the nav team's NavigationDictionary.sldd. The
% harness links our own dictionary (config/hopper_harness.sldd) that
% references theirs, with base-workspace access on, so buildSimInput's
% variables still reach the model and their file is never edited.
dd = fullfile(hopperRoot('config'), 'hopper_harness.sldd');
if ~isfile(dd)
    if ~isfolder(fileparts(dd)), mkdir(fileparts(dd)); end
    d = Simulink.data.dictionary.create(dd);
    addDataSource(d, 'NavigationDictionary.sldd');
    d.EnableAccessToBaseWorkspace = true;
    saveChanges(d);
    close(d);
end
set_param(m, 'DataDictionary', 'hopper_harness.sldd');
end

function setScript(blk, lines, params, sampleTime)
% MATLAB Function block code, workspace parameters and sample time
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
ch.Script = strjoin(lines, newline);
for p = params
    d = ch.find('-isa', 'Stateflow.Data', 'Name', p{1});
    d.Scope = 'Parameter';
    try, d.Tunable = false; catch, end   % structs as fixed parameters
end
set_param(blk, 'SystemSampleTime', sampleTime);   % the chart's own SampleTime property doesn't stick
if ~strcmp(get_param(blk, 'SystemSampleTime'), sampleTime)
    error('build_hopper_6dof:sampleTime', '%s: sample time %s did not take.', blk, sampleTime);
end
end

function h = port(m, blk, k)
ph = get_param([m '/' blk], 'PortHandles');
h = ph.Outport(k);
end

function s = legacyGainSelector()
% Gain_Selector from archive/models/hopper_6dof_lqr.slx and LQ_6DOF_V3.slx
% (March): v2's scheduler without the quaternion sign fix
s = strjoin({
    'function u = Gain_Selector(t, x, tgrid, xref, K1flat, K2grid, unom)'
    '%#codegen'
    '% Legacy (March) scheduler: no quaternion hemisphere check'
    'u = zeros(4,1);'
    'x = x(:);'
    'N = length(tgrid);'
    'if t <= tgrid(1)'
    '    idx = 1;'
    '    w = 0.0;'
    'elseif t >= tgrid(end)'
    '    idx = N-1;'
    '    w = 1.0;'
    'else'
    '    idx = 1;'
    '    for k = 1:(N-1)'
    '        if tgrid(k) <= t && t < tgrid(k+1)'
    '            idx = k;'
    '            break;'
    '        end'
    '    end'
    '    t1 = tgrid(idx);'
    '    t2 = tgrid(idx+1);'
    '    w  = (t - t1) / (t2 - t1);'
    'end'
    'xref_t = (1-w)*xref(idx,:)'' + w*xref(idx+1,:)'';'
    'Kflat = (1-w)*K1flat(idx,:) + w*K1flat(idx+1,:);'
    'K1 = reshape(Kflat, 13, 4)'';'
    'K2 = (1-w)*K2grid(idx,:)'' + w*K2grid(idx+1,:)'';'
    'u_nom = (1-w)*unom(idx,:)'' + w*unom(idx+1,:)'';'
    'dx = (x - xref_t);'
    'u  = u_nom - K1*dx - K2;'
    'end'
    }, newline);
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

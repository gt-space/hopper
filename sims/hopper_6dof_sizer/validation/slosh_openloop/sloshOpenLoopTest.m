function R = sloshOpenLoopTest(slxFile, ox_frac, fu_frac, inputMode, floorSet)
% Open-loop free-oscillation / sign test of the Slosh subsystem.
% The Slosh subsystem of slxFile drives a planar rigid body (lateral y,
% yaw about body z) with constant thrust acceleration f along body x
% (hover-equivalent). A 50 N, 50 ms lateral pulse at the CG excites it.
% Truth: independent Newton-Euler model of the dry vehicle + pendulum bob
% (SP-106 analog) on a frictionless hinge. Damping set to zero.
%   inputMode 'specific'  : pendulum driven by specific force (fixed wiring)
%             'kinematic' : driven by dVb/dt incl. gravity (old wiring)
%   floorSet  'asbuilt' | 'sizer' : tank floor used by the truth model
% Run sim_setup first (uses OUT, MoI_init, cg_init, TANKS from the base
% workspace). See SLOSH_VALIDATION.md section 7.
%
%   R = sloshOpenLoopTest('hopper_6dof_NED_v2.slx', 1, 0, 'specific', 'asbuilt')
for v = {'OUT', 'MoI_init', 'cg_init', 'TANKS'}, W.(v{1}) = evalin('base', v{1}); end
[~, src] = fileparts(slxFile);
load_system(slxFile);

g = 9.81; f = g; R_t = W.TANKS.singular.radius; xi = 1.841;
Mv = W.OUT.Vehicle.WetMass; Izz = W.MoI_init(3,3); cg = W.cg_init;
F0 = 50; tp = 0.05; Tend = 6; dt = 1e-4;

assignin('base', 'TANKS', W.TANKS);
assignin('base', 'ox_mass_profile', linspace(0.1, 16, 200));
assignin('base', 'fu_mass_profile', linspace(0.1, 13, 200));
for v = {'ox_lateral_damping_ratios', 'ox_axial_damping_ratios', 'fu_lateral_damping_ratios', 'fu_axial_damping_ratios'}
    assignin('base', v{1}, zeros(1, 200));
end
m_ox = ox_frac * 16; m_fu = fu_frac * 13;   % unused tank: 0 (saturation lowered below)

h = 'sloshHarness';
if bdIsLoaded(h), close_system(h, 0); end
new_system(h);
add_block([src '/Slosh'], [h '/Slosh'], 'Position', [400 100 560 300]);
set_param([h '/Slosh/Saturation'], 'LowerLimit', '0'); set_param([h '/Slosh/Saturation1'], 'LowerLimit', '0');
add_block('simulink/Sources/Constant', [h '/mox'], 'Value', num2str(m_ox, 17), 'Position', [250 100 280 120]);
add_block('simulink/Sources/Constant', [h '/mfu'], 'Value', num2str(m_fu, 17), 'Position', [250 130 280 150]);
add_block('simulink/Sources/Constant', [h '/cg'], 'Value', sprintf('[%.17g 0 0]', cg), 'Position', [250 280 300 300]);
add_block('simulink/Sources/Constant', [h '/P'], 'Value', ...
    sprintf('[%.17g %.17g %.17g %.17g %.17g %d %.17g]', Mv, Izz, f, F0, tp, strcmp(inputMode, 'kinematic'), g), ...
    'Position', [50 400 200 420]);
add_block('simulink/Sources/Clock', [h '/t'], 'Position', [50 450 70 470]);
add_block('simulink/User-Defined Functions/MATLAB Function', [h '/Plant'], 'Position', [250 380 350 520]);
rt = sfroot;
ch = rt.find('-isa', 'Stateflow.EMChart', 'Path', [h '/Plant']);
ch.Script = sprintf([ ...
 'function [xdot, a_in, w, wdot] = fcn(x, F, Mo, t, P)\n' ...
 'M = P(1); Izz = P(2); f = P(3); Fext = P(4) * (t < P(5));\n' ...
 'xdot = [x(2); (M*f*x(3) + F(2) + Fext)/M; x(4); Mo(3)/Izz];\n' ...
 'a_in = [f; (F(2) + Fext)/M; F(3)/M];\n' ...
 'if P(6) > 0  %% old wiring: F_total/m incl. gravity (hover thrust = weight)\n' ...
 '    a_in = a_in + [-P(7); P(7)*x(3); 0];\n' ...
 'end\n' ...
 'w = [0; 0; x(4)]; wdot = [0; 0; Mo(3)/Izz];\n']);
add_block('simulink/Continuous/Integrator', [h '/x'], 'InitialCondition', 'zeros(4,1)', 'Position', [400 400 430 430]);
add_block('simulink/Discrete/Memory', [h '/memF'], 'InitialCondition', 'zeros(3,1)', 'Position', [620 100 650 130]);
add_block('simulink/Discrete/Memory', [h '/memM'], 'InitialCondition', 'zeros(3,1)', 'Position', [620 150 650 180]);
tw = {'X', 'x/1'; 'F', 'Slosh/1'; 'Mo', 'Slosh/2'; 'Pox', 'Slosh/3'; 'Pfu', 'Slosh/5'};
for i = 1:size(tw, 1)
    add_block('simulink/Sinks/To Workspace', [h '/tw' tw{i,1}], 'VariableName', ['tw' tw{i,1}], ...
        'SaveFormat', 'Array', 'Position', [750 50+50*i 800 70+50*i]);
    add_line(h, tw{i,2}, ['tw' tw{i,1} '/1']);
end
for k = [4 6]
    add_block('simulink/Sinks/Terminator', sprintf('%s/term%d', h, k), 'Position', [620 200+20*k 640 210+20*k]);
    add_line(h, sprintf('Slosh/%d', k), sprintf('term%d/1', k));
end
add_line(h, 'mox/1', 'Slosh/1'); add_line(h, 'mfu/1', 'Slosh/2'); add_line(h, 'cg/1', 'Slosh/6');
add_line(h, 'Plant/2', 'Slosh/3'); add_line(h, 'Plant/3', 'Slosh/4'); add_line(h, 'Plant/4', 'Slosh/5');
add_line(h, 'Slosh/1', 'memF/1'); add_line(h, 'Slosh/2', 'memM/1');
add_line(h, 'x/1', 'Plant/1'); add_line(h, 'memF/1', 'Plant/2'); add_line(h, 'memM/1', 'Plant/3');
add_line(h, 't/1', 'Plant/4'); add_line(h, 'P/1', 'Plant/5'); add_line(h, 'Plant/1', 'x/1');
set_param(h, 'SolverType', 'Fixed-step', 'Solver', 'ode4', 'FixedStep', num2str(dt), 'StopTime', num2str(Tend), ...
    'ReturnWorkspaceOutputs', 'on', 'AlgebraicLoopMsg', 'error');
so = sim(h);
for v = {'twX','twF','twMo','twPox','twPfu'}, A = so.(v{1}); if ndims(A) == 3, A = squeeze(A).'; end; sig.(v{1}) = A; end
t = (0:dt:Tend)'; n = min(numel(t), size(sig.twX, 1)); t = t(1:n);
Y = sig.twX(1:n, 1); th = sig.twX(1:n, 3);
if ox_frac >= fu_frac, ysl = sig.twPox(1:n, 2); else, ysl = sig.twPfu(1:n, 2); end

%% Truth: dry body + one pendulum (the excited tank), Newton-Euler, linear
if ox_frac >= fu_frac, rho = 1141; m = m_ox; Hc = 29.18; Hs = W.TANKS.singular.oxidizer.h; zc = 1.854;
else,                  rho = 800;  m = m_fu; Hc = 18.91; Hs = W.TANKS.singular.fuel.h;     zc = 1.0668; end
if strcmp(floorSet, 'asbuilt'), floor_ = zc - Hc*0.0254/2; else, floor_ = zc - Hs/2; end
hl = m / rho / (pi * R_t^2);
m1 = m * 2 * R_t * tanh(xi*hl/R_t) / (xi * (xi^2 - 1) * hl);
L  = R_t / (xi * tanh(xi*hl/R_t));
H1 = hl - (2*R_t/xi) * tanh(xi*hl/(2*R_t));
xH1 = floor_ + H1 - cg; xh = xH1 + L;
MD = Mv - m1; rD = -m1 * xH1 / MD; ID = Izz - m1*xH1^2 - MD*rD^2;
% q = [Y; th; psi]; Mq'' = K q + B Fext
Mm = [m1, m1*xH1, -m1*L; MD, MD*rD, 0; 0, ID, 0];
K  = [0, m1*f, m1*f; 0, MD*f, -m1*f; 0, 0, -(xh - rD)*m1*f];
B  = [0; 1; -rD];
rhs = @(tt, z) [z(4:6); Mm \ (K*z(1:3) + B*F0*(tt < tp))];
[~, Z] = ode45(rhs, t, zeros(6,1), odeset('RelTol', 1e-10, 'AbsTol', 1e-12, 'MaxStep', 1e-3));
Yt = Z(:,1); tht = Z(:,2); yt = -L*Z(:,3);

w_unc = sqrt(f/L);
R = struct('src', src, 'mode', inputMode, 'ox', ox_frac, 'fu', fu_frac, 'm1', m1, 'L', L, 'xH1', xH1, ...
    'w_uncoupled', w_unc, 'w_truth', domFreq(t, yt), 'w_model', domFreq(t, ysl), ...
    'th_peak_truth', max(abs(tht)), 'th_peak_model', max(abs(th)), ...
    'y_peak_truth', max(abs(yt)), 'y_peak_model', max(abs(ysl)), ...
    'err_th', norm(th - tht) / norm(tht), 'err_y', norm(ysl - yt) / norm(yt), 'err_Y', norm(Y - Yt) / norm(Yt), ...
    'F0sign', NaN);
% Sign check: first sample of the pulse, pendulum force on the vehicle
% should be +m1*a (liquid does not follow yet), i.e. along the pulse.
Fy = sig.twF(1:n, 2); i1 = find(t > 0.002, 1);
R.F0sign = Fy(i1) / (m1 * F0 / Mv);
R.t = t; R.th = th; R.tht = tht; R.y = ysl; R.yt = yt; R.Y = Y; R.Yt = Yt;
close_system(h, 0);
end

function w = domFreq(t, y)
y = y(t > 0.2); t = t(t > 0.2); y = y - mean(y);
N = 2^nextpow2(16*numel(y)); Yf = abs(fft(y .* hann(numel(y)), N));
fr = (0:N-1)' / N / (t(2) - t(1)) * 2*pi;
sel = fr > 1 & fr < 100; [~, k] = max(Yf .* sel); w = fr(k);
end

function verify_env_c()
%VERIFY_ENV_C  Check the generated hopper_env C code against Simulink.
%   Builds a gust-free copy of the environment C code (build/verify_nogust),
%   runs hopper_6dof_NED_v2_fswBridge in Simulink with gusts off, records the
%   controller's commands and the Environment outputs, replays the same
%   commands through the C code (replay_main.c), and prints the differences.
%
%   Gusts are off on both sides because the generated code draws different
%   random numbers than Simulink's interpreted gust blocks. Run in a fresh
%   MATLAB session (build_env_c refuses if hopper_env is already loaded).
%   The bridge model is never saved.

here  = fileparts(mfilename('fullpath'));
sizer = fileparts(here);
build = fullfile(here, 'build', 'verify_nogust');
cd(sizer);

if ~evalin('base', 'exist(''IN'', ''var'') && exist(''VEH'', ''var'')')
    evalin('base', 'sim_setup_cached');
end
code = build_env_c('Gusts', false, 'Folder', build);

%% 1. Simulink run (gusts off) with the command and x_true signals logged
mdl = 'hopper_6dof_NED_v2_fswBridge';
env = [mdl '/Environment'];
load_system(mdl);
set_param([env '/Bernoulli Binary Generator'], 'ProbabilityOfZero', '1');
cph = get_param([mdl '/Controller'], 'PortHandles');
set_param(cph.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'u_cmd');
eph = get_param(env, 'PortHandles');
xPort = str2double(get_param([env '/x_true'], 'Port'));
set_param(eph.Outport(xPort), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'x_true');
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'logsout');

fprintf('Running Simulink reference...\n');
out = sim(Simulink.SimulationInput(mdl));
u  = out.logsout.get('u_cmd').Values;
xs = out.logsout.get('x_true').Values;
t  = u.Time;
U  = asRows(u.Data, numel(t));
Xs = asRows(xs.Data, numel(xs.Time));
thrustSL = asRows(out.get('thrust').Data, numel(out.get('thrust').Time));
zSL      = asRows(out.get('z').Data, numel(out.get('z').Time));

cmdFile = fullfile(build, 'replay_commands.csv');
fid = fopen(cmdFile, 'w');
fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g\n', [t U].');
fclose(fid);

%% 2. Link the replay program against the already-compiled model objects
cc  = mex.getCompilerConfigurations('C', 'Selected');
bin = fullfile(cc.Location, 'bin');
gcc = fullfile(bin, 'gcc.exe'); gxx = fullfile(bin, 'g++.exe');
inc = sprintf('-I"%s" -I"%s" -I"%s" -I"%s" -I"%s"', code, build, ...
    fullfile(matlabroot, 'extern', 'include'), fullfile(matlabroot, 'simulink', 'include'), ...
    fullfile(matlabroot, 'rtw', 'c', 'src'));
defs = ['-DCLASSIC_INTERFACE=0 -DALLOCATIONFCN=0 -DTERMFCN=1 -DONESTEPFCN=1 -DMAT_FILE=0 ' ...
        '-DMULTI_INSTANCE_CODE=0 -DINTEGER_CODE=0 -DMT=0 -DTID01EQ=1 -DMODEL=hopper_env -DHAVESTDIO'];
replayObj = fullfile(build, 'replay_main.obj');
exe       = fullfile(build, 'replay.exe');
objs = dir(fullfile(code, '*.obj'));
objs = objs(~strcmp({objs.name}, 'ert_main.obj'));
objList = strjoin(cellfun(@(n) ['"' fullfile(code, n) '"'], {objs.name}, 'UniformOutput', false), ' ');

run_or_fail(sprintf('"%s" -c -fwrapv -m64 -O0 -msse2 %s %s -o "%s" "%s"', gcc, defs, inc, replayObj, fullfile(here, 'replay_main.c')));
run_or_fail(sprintf('"%s" -static -m64 -o "%s" "%s" %s -lws2_32', gxx, exe, replayObj, objList));

%% 3. Replay through C and compare
cOut = fullfile(build, 'replay_outputs.csv');
run_or_fail(sprintf('"%s" "%s" "%s"', exe, cmdFile, cOut));
C = readmatrix(cOut);

% Each C row is labelled with the model time after its step; pair it with
% the Simulink sample at that same time.
dt = t(2) - t(1);
k  = round(C(:, 1) / dt) + 1;
keep = k <= min([size(Xs,1), size(thrustSL,1), size(zSL,1)]);
C = C(keep, :); k = k(keep);
S = [Xs(k, :), thrustSL(k), zSL(k)];
tc = C(:, 1);

names = [{'pos_n','pos_e','pos_d','vel_n','vel_e','vel_d','p','q','r','q0','q1','q2','q3'}, {'thrust','z'}];
fprintf('\nC vs Simulink over %d steps (t = %.3f .. %.3f s):\n', numel(tc), tc(1), tc(end));
fprintf('%-8s %14s %14s\n', 'signal', 'max|diff|', 'max|value|');
for i = 1:numel(names)
    fprintf('%-8s %14.3g %14.4g\n', names{i}, max(abs(C(:, 1+i) - S(:, i))), max(abs(S(:, i))));
end

% Replay is open loop, so an unstable plant amplifies any tiny difference.
% When differences first appear says more than their final size.
D = max(abs(C(:, 2:14) - S(:, 1:13)), [], 2);
fprintf('\nFirst time the largest x_true difference exceeds:\n');
for tol = [1e-12 1e-9 1e-6 1e-3 1]
    j = find(D > tol, 1);
    if isempty(j), fprintf('  %-6g never\n', tol); else, fprintf('  %-6g t = %.3f s\n', tol, tc(j)); end
end

% Sensors: compare at each sensor's own logged sample times, up to the
% point the state has drifted (after that the inputs differ anyway).
tDrift = tc(find([D; inf] > 1e-9, 1));
sensors = {'imu', 6; 'mag', 3; 'baro', 2; 'gps', 8; 'lidar', 4};
col = 17;
fprintf('\nSensors, C vs Simulink, up to t = %.3f s:\n', tDrift);
for s = 1:size(sensors, 1)
    ts = out.get([sensors{s, 1} '_meas']);
    V  = asRows(ts.Data, numel(ts.Time));
    j  = round(ts.Time / dt) + 1;                 % Simulink sample -> C row
    [ok, r] = ismember(j, k);
    ok = ok & ts.Time <= tDrift;
    d = max(abs(C(r(ok), col:col + sensors{s, 2} - 1) - V(ok, :)), [], 'all');
    fprintf('  %-6s %6d samples  max|diff| %.3g\n', sensors{s, 1}, nnz(ok), d);
    col = col + sensors{s, 2};
end
save(fullfile(build, 'verify_env_c.mat'), 't', 'U', 'Xs', 'C', 'thrustSL', 'zSL');
end

function A = asRows(D, nT)
% Logged data as one row per time step
if ndims(D) == 3
    A = reshape(permute(D, [3 1 2]), nT, []);
elseif size(D, 1) == nT
    A = D;
else
    A = D.';
end
end

function run_or_fail(cmd)
[status, msg] = system(cmd);
if ~isempty(strtrim(msg)), fprintf('%s\n', strtrim(msg)); end
if status ~= 0
    error('verify_env_c:cmd', 'Command failed (%d): %s', status, cmd);
end
end

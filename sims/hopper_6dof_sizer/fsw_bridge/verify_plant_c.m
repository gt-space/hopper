function R = verify_plant_c(varargin)
%VERIFY_PLANT_C  Check the generated plant against Simulink, open loop.
%   R = verify_plant_c() flies the CTV in Simulink with its own LQR
%   (CONTROLLER=lqr_sched), samples the command at the flight software's
%   period, and replays that command history twice:
%     1. in Simulink with CONTROLLER=fsw (the same harness, command on the
%        root inport u_fsw)
%     2. in the generated C (build_plant_c, x86_64, gusts off), compiled
%        and run in WSL with plant_replay.c
%   and compares x_fsw at every command. Gusts are off in all three runs,
%   since generated code draws different random numbers than Simulink.
%   It also compares the sensor ports (sens_meas, sens_valid, sens_stamp)
%   at every command.
%
%   Options: 'Vehicle' (default 'ctv'), 'Build' (default true: regenerate
%   the code; false reuses build/verify/plant), 'Tol' (default 1e-6, max
%   abs difference allowed in any state), 'Excite' (default true: add
%   0.5 deg TVC and 2 mN m roll sines to the replayed command, so the
%   lateral and attitude states move too), 'SensTol' (default 1e-6, max
%   abs difference allowed in any sensor channel, valid flag or stamp).
%
%   R has the max abs error per state, the time of the worst one, both
%   state histories, the sensor comparison (R.sens) and pass. Needs WSL with gcc (wsl -e gcc).

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Vehicle', 'ctv');
p.addParameter('Build', true, @islogical);
p.addParameter('Tol', 1e-6, @isnumeric);
p.addParameter('Excite', true, @islogical);
p.addParameter('SensTol', 1e-6, @isnumeric);
p.parse(varargin{:});
vehicle = char(p.Results.Vehicle);
dir_ = fullfile(here, 'build', 'verify');
if ~isfolder(dir_), mkdir(dir_); end
m = 'hopper_6dof';

% ---- 1. reference command: the harness's own LQR, gusts off ----
logFile = fullfile(dir_, 'reference_log.mat');
Vlqr = variantConfig('VEHICLE', vehicle, 'LOG', logFile);   % the signal dump has gnc_bus.u_cmd
[in, W] = buildSimInput(m, [], [], Vlqr);
in = noGusts(in, m);
ref = sim(in);
dt = 0.01;   % the flight software's period (fswWorkspace)
L = load(logFile, 'logsout');
ucmd = L.logsout.get('gnc_bus').Values.u_cmd;
u = squeeze(ucmd.Data);
if size(u, 1) == 4 && size(u, 2) ~= 4, u = u.'; end
tc = (0:dt:ref.tout(end) - dt).';
uc = interp1(ucmd.Time, u, tc, 'previous', 'extrap');
if p.Results.Excite
    % calm air and no gusts fly straight up, which leaves the lateral and
    % attitude states at zero; small gimbal and roll commands exercise them
    uc(:, 2) = uc(:, 2) + deg2rad(0.5) * sin(2*pi*0.5*tc);
    uc(:, 3) = uc(:, 3) + deg2rad(0.5) * cos(2*pi*0.3*tc);
    uc(:, 4) = uc(:, 4) + 2e-3 * sin(2*pi*0.2*tc);
end
cmd = [tc, uc];
fid = fopen(fullfile(dir_, 'commands.csv'), 'w');
fprintf(fid, '%.17g,%.17g,%.17g,%.17g,%.17g\n', cmd.');
fclose(fid);
fprintf('Reference LQR flight: %.2f s, %d commands at %g s\n', ref.tout(end), numel(tc), dt);

% ---- 2. the same commands through CONTROLLER=fsw in Simulink ----
Vfsw = variantConfig('VEHICLE', vehicle, 'CONTROLLER', 'fsw');
W.fsw_u = cmd;
in = noGusts(buildSimInput(m, [], W, Vfsw), m);
sl = sim(in);
xs = rootPort(sl, 'x_fsw');
xsl = squeeze(xs.Data);
if size(xsl, 1) == 13 && size(xsl, 2) ~= 13, xsl = xsl.'; end
xsl = interp1(xs.Time, xsl, tc, 'previous');
ssl = [portAt(sl, 'sens_meas', 26, tc), portAt(sl, 'sens_valid', 9, tc), portAt(sl, 'sens_stamp', 9, tc)];

% ---- 3. the generated C ----
code = fullfile(dir_, 'plant');
if p.Results.Build
    close_system(m, 0);   % build_plant_c loads its own copy
    build_plant_c('Variants', Vfsw, 'Gusts', false, 'Target', 'x86_64', 'Folder', code);
end
info = jsondecode(fileread(fullfile(code, 'plant.json')));
srcs = readlines(fullfile(code, 'sources.txt'));
defs = readlines(fullfile(code, 'defines.txt'));
srcs = srcs(strlength(srcs) > 0); defs = defs(strlength(defs) > 0);
w = @(f) wslPath(f);
cc = sprintf('gcc -O2 -o %s %s %s -I%s %s -lm', w(fullfile(dir_, 'plant_replay')), ...
    w(fullfile(here, 'plant_replay.c')), strjoin(w(fullfile(code, 'src', srcs)), ' '), ...
    w(fullfile(code, 'src')), strjoin("-D" + defs, ' '));
run(sprintf('wsl -e bash -c "%s"', cc));
run(sprintf('wsl -e %s %s %s %d %s', w(fullfile(dir_, 'plant_replay')), w(fullfile(dir_, 'commands.csv')), ...
    w(fullfile(dir_, 'states_c.csv')), round(dt / info.fixed_step), w(fullfile(dir_, 'sensors_c.csv'))));
xc = readmatrix(fullfile(dir_, 'states_c.csv'));
tcC = xc(:, 1); xc = xc(:, 2:14);
sc = readmatrix(fullfile(dir_, 'sensors_c.csv'));
sc = sc(:, 2:end);

% ---- compare where both ran ----
n = min([size(xc, 1), find(~isnan(xsl(:, 1)), 1, 'last')]);
err = abs(xc(1:n, :) - xsl(1:n, :));
[worst, at] = max(err, [], 1);
serr = abs(sc(1:n, :) - ssl(1:n, :));
serr(isnan(sc(1:n, :)) & isnan(ssl(1:n, :))) = 0;
[sworst, sat] = max(serr, [], 1);
L = sensor_layout();
groups = {'accel', L.ACCEL; 'gyro', L.GYRO; 'baro', L.BARO; 'lidar', L.LIDAR; 'gps_pos', L.GPS_POS; ...
    'gps_fix/sats/age', [L.GPS_FIX L.GPS_SATS L.GPS_AGE]; 'mag', L.MAG; 'gps_vel', L.GPS_VEL; ...
    'dgps', L.DGPS; 'valid', 26 + (1:9); 'stamp', 35 + (1:9)};
sens = struct('max_err', sworst, 'worst_time', tcC(sat).', 's_c', sc(1:n, :), 's_sim', ssl(1:n, :), ...
    'pass', all(sworst < p.Results.SensTol));
R = struct('max_err', worst, 'worst_time', tcC(at).', 'n', n, 'x_c', xc(1:n, :), ...
    'x_sim', xsl(1:n, :), 't', tcC(1:n), 'sens', sens, ...
    'pass', all(worst < p.Results.Tol) && sens.pass);
names = {'pN', 'pE', 'pD', 'vN', 'vE', 'vD', 'p', 'q', 'r', 'q0', 'q1', 'q2', 'q3'};
fprintf('C vs Simulink over %d commands (%.2f s):\n', n, tcC(n));
for k = 1:13, fprintf('  %-3s max |err| %.3g at t = %.2f\n', names{k}, worst(k), R.worst_time(k)); end
fprintf('Sensor ports:\n');
for g = 1:size(groups, 1)
    fprintf('  %-16s max |err| %.3g\n', groups{g, 1}, max(sworst(groups{g, 2})));
end
fprintf('  (sensor tol %g)\n', p.Results.SensTol);
fprintf('%s (tol %g)\n', ternary(R.pass, 'PASS', 'FAIL'), p.Results.Tol);
end

function v = portAt(out, name, n, t)
% a root outport's value at times t (held), n columns
s = rootPort(out, name);
d = squeeze(double(s.Data));
if size(d, 1) == n && size(d, 2) ~= n, d = d.'; end
v = interp1(s.Time, d, t, 'previous');
end

function s = rootPort(out, name)
% a root outport's timeseries from yout (its elements are unnamed: match the block)
for k = 1:out.yout.numElements
    e = out.yout{k};
    if strcmp(e.BlockPath.getBlock(1), ['hopper_6dof/' name]), s = e.Values; return; end
end
error('verify_plant_c:port', 'no root outport %s in yout', name);
end

function in = noGusts(in, m)
load_system(m);
gate = find_system(m, 'MatchFilter', @Simulink.match.allVariants, 'LookUnderMasks', 'all', ...
    'FollowLinks', 'on', 'Name', 'Bernoulli Binary Generator');
for k = 1:numel(gate), in = in.setBlockParameter(gate{k}, 'ProbabilityOfZero', '1'); end
end

function s = wslPath(f)
f = string(f);
s = "/mnt/" + lower(extractBefore(f, 2)) + replace(extractAfter(f, 2), "\", "/");
end

function run(cmd)
[status, outText] = system(cmd);
fprintf('%s', outText);
if status ~= 0, error('verify_plant_c:cmd', 'failed (%d): %s', status, cmd); end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function outDir = export_controller_vectors(varargin)
%EXPORT_CONTROLLER_VECTORS  Golden (t, x) -> u vectors for the flight software's LQR.
%   export_controller_vectors() writes, next to the gain tables
%   export_plant_gains writes (fsw_bridge/build/gains/<vehicle>):
%     controller_vectors.csv  header t,x1..x13,u1..u4; one case per row,
%                             full double precision
%     controller_vectors.json what the cases are, the tolerance, and the
%                             SHA-256 of the gain CSVs they were made with
%   The flight software's controller, loaded with those CSVs, must
%   reproduce every u to within the tolerance (1e-9 relative to max(1, |u|)).
%
%   u comes from the model's own scheduler code: the Gain_Selector MATLAB
%   Function in the vehicle controller's subsystem file
%   (subsystems/gnc/LQR_Gain_Scheduler.slx, which both lqr_sched and
%   LQR_CTV_Tuning.slx reference) is read out of the file and run
%   as is, with the exported tables. So the vectors test the law, the
%   interpolation, the clamping, the K1 layout and the quaternion
%   hemisphere flip exactly as Simulink flies them.
%
%   Cases:
%     flight    the nominal flight's state every 10 ms (the flight
%               software's period), at its own time
%     grid      every 25th tgrid point, exactly on the grid, x = xref
%     dispersed random states around the reference (seed 1): 0.5 m,
%               0.5 m/s, 0.5 rad/s, 0.1 in each quaternion component
%               (renormalized), half of them with the quaternion negated
%     ends      times before the start and after the end of the schedule
%
%   Options: 'Vehicle' (default 'ctv'), 'Controller' (default the
%   vehicle's), 'Folder' (default fsw_bridge/build/gains/<vehicle>).

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Vehicle', 'ctv');
p.addParameter('Controller', '');
p.addParameter('Folder', '');
p.parse(varargin{:});
vehicle = char(p.Results.Vehicle);
outDir = export_plant_gains('Vehicle', vehicle, 'Controller', p.Results.Controller, ...
    'Folder', p.Results.Folder);
info = jsondecode(fileread(fullfile(outDir, 'gains.json')));

% ---- the tables exactly as exported ----
rd = @(n) readmatrix(fullfile(outDir, [n '.csv']));
tgrid = rd('tgrid'); xref = rd('xref'); k1flat = rd('k1flat'); k2grid = rd('k2grid'); unom = rd('unom');

% ---- the scheduler's code, out of the subsystem file ----
switch info.controller
    case "lqr_ctv_tuning", file = 'LQR_CTV_Tuning';
    case "lqr_sched",      file = 'LQR_Gain_Scheduler';
end
[gainSelector, srcSha, where] = schedulerFunction(file);

% ---- cases ----
% the controller's own input, nav_bus.x_nav, from the signal dump
logFile = fullfile(tempdir, 'hopper_controller_vectors', 'flight_log.mat');
V = variantConfig('VEHICLE', vehicle, 'CONTROLLER', char(info.controller), 'LOG', logFile);
sim(buildSimInput('hopper_6dof', [], [], V));
L = load(logFile, 'logsout');
xs = L.logsout.get('nav_bus').Values.x_nav;
X = squeeze(xs.Data);
if size(X, 2) ~= 13, X = X.'; end
assert(size(X, 2) == 13, 'export_controller_vectors:x', 'logged x is not the 13-state');
tf = (0:0.01:xs.Time(end)).';
Xf = interp1(xs.Time, X, tf, 'previous');
cases = {'flight', tf, Xf};

ig = (1:25:numel(tgrid)).';
cases(end+1, :) = {'grid', tgrid(ig), xref(ig, :)};

rng(1);
nd = 400;
td = tgrid(1) + (tgrid(end) - tgrid(1)) * rand(nd, 1);
Xd = interp1(tgrid, xref, td);
Xd(:, 1:3) = Xd(:, 1:3) + 0.5 * randn(nd, 3);
Xd(:, 4:6) = Xd(:, 4:6) + 0.5 * randn(nd, 3);
Xd(:, 7:9) = Xd(:, 7:9) + 0.5 * randn(nd, 3);
q = Xd(:, 10:13) + 0.1 * randn(nd, 4);
q = q ./ vecnorm(q, 2, 2);
flip = (1:nd).' > nd / 2;
q(flip, :) = -q(flip, :);
Xd(:, 10:13) = q;
cases(end+1, :) = {'dispersed', td, Xd};

te = [tgrid(1) - 1; tgrid(1) - 0.01; tgrid(end) + 0.01; tgrid(end) + 5];
cases(end+1, :) = {'ends', te, xref([1 1 end end], :) + 0.01};

% ---- u from the model's code ----
rows = [];
counts = struct();
for c = 1:size(cases, 1)
    t = cases{c, 2}; Xc = cases{c, 3};
    U = zeros(numel(t), 4);
    for k = 1:numel(t)
        U(k, :) = gainSelector(t(k), Xc(k, :).', tgrid, xref, k1flat, k2grid, unom).';
    end
    rows = [rows; t, Xc, U]; %#ok<AGROW>
    counts.(cases{c, 1}) = numel(t);
end
assert(all(isfinite(rows(:))), 'export_controller_vectors:finite', 'non-finite vector');

f = fullfile(outDir, 'controller_vectors.csv');
fid = fopen(f, 'w');
fprintf(fid, '%s\n', strjoin([{'t'}, compose('x%d', 1:13), compose('u%d', 1:4)], ','));
fprintf(fid, [strjoin(repmat({'%.17g'}, 1, 18), ',') '\n'], rows.');
fclose(fid);

vinfo = struct('vehicle', vehicle, 'controller', info.controller, ...
    'source', where, ...
    'source_sha256', srcSha, 'gains_sha256', info.sha256, 'git_commit', info.git_commit, ...
    'columns', 't s; x1..x13 pos NED m, vel NED m/s, body rates rad/s, quaternion scalar-first; u1..u4 thrust N, TVC pitch rad, TVC yaw rad, RCS torque N m', ...
    'tolerance', '|u - u_ref| <= 1e-9 * max(1, |u_ref|) per component', ...
    'cases', counts, 'rows', size(rows, 1), ...
    'exported', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
fid = fopen(fullfile(outDir, 'controller_vectors.json'), 'w');
fprintf(fid, '%s\n', jsonencode(vinfo, 'PrettyPrint', true));
fclose(fid);
fprintf('%d controller vectors (%s) -> %s\n', size(rows, 1), info.controller, f);
end

function [fh, sha, where] = schedulerFunction(file)
% The Gain_Selector MATLAB Function's code out of a subsystem file (or a
% subsystem file it references, e.g. LQR_CTV_Tuning -> LQR_Gain_Scheduler),
% as a callable function handle, and the SHA-256 of that code.
[code, where] = findGainSelector(file, {});
assert(~isempty(code), 'export_controller_vectors:source', 'no Gain_Selector in %s', file);
fprintf('Gain_Selector from %s\n', where);
d = fullfile(tempdir, 'hopper_controller_vectors');
if ~isfolder(d), mkdir(d); end
fid = fopen(fullfile(d, 'Gain_Selector.m'), 'w');
fprintf(fid, '%s', code);
fclose(fid);
addpath(d);
clear Gain_Selector
fh = @Gain_Selector;
sha = sha256File(fullfile(d, 'Gain_Selector.m'));
end

function [code, where] = findGainSelector(file, seen)
code = ''; where = '';
if ismember(file, seen), return; end
seen{end+1} = file;
wasLoaded = bdIsLoaded(file);
load_system(file);
if ~wasLoaded, closeIt = onCleanup(@() close_system(file, 0)); end
charts = sfroot().find('-isa', 'Stateflow.EMChart');
for c = charts(:)'
    if startsWith(c.Path, [file '/']) && contains(c.Script, 'function u = Gain_Selector')
        code = c.Script; where = c.Path; return
    end
end
subs = find_system(file, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'BlockType', 'SubSystem');
for k = 1:numel(subs)
    ref = get_param(subs{k}, 'ReferencedSubsystem');
    if ~isempty(ref)
        [code, where] = findGainSelector(ref, seen);
        if ~isempty(code), return; end
    end
end
end

function outDir = export_plant_gains(varargin)
%EXPORT_PLANT_GAINS  Write a vehicle's LQR gain schedule for the flight software.
%   export_plant_gains() writes the time-scheduled LQR tables the CTV flies
%   in models/hopper_6dof.slx to fsw_bridge/build/gains/ctv as CSV (one row
%   per grid time, full double precision, no header) plus gains.json saying
%   what they are:
%     tgrid.csv   n x 1    s, schedule times from the start of the flight
%     xref.csv    n x 13   reference state (Target_Trajectory): pos NED m,
%                          vel NED m/s, body rates rad/s, quaternion
%                          scalar-first
%     k1flat.csv  n x 52   K1 (4 x 13) per row, ROW-major: the 13 thrust
%                          gains, then TVC pitch, TVC yaw, RCS
%                          (flat(13*i + j) = K1(i, j), 0-based), which is
%                          MATLAB's reshape(K1.', 1, 52)
%     k2grid.csv  n x 4    K2 per row
%     unom.csv    n x 4    nominal command per row: thrust N, TVC pitch
%                          rad, TVC yaw rad, RCS torque N m
%   The harness flies u = unom - K1 * (x - xref) - K2, every table
%   linearly interpolated in t and held at its ends, with x's quaternion
%   negated when dot(q, qref) < 0 (Gain_Selector in the LQR subsystems).
%
%   The tables are the ones the vehicle's controller reads: the CTV
%   defaults to CONTROLLER=lqr_ctv_tuning, whose gains are K1flat_ctv and
%   K2grid_ctv (src/gnc/ctvTuningGains.m); lqr_sched reads K1flat and
%   K2grid. Both share tgrid, xref and unom.
%
%   Options: 'Vehicle' (default 'ctv'), 'Controller' (default the
%   vehicle's default, lqr_ctv_tuning for the CTV; lqr_sched or
%   lqr_ctv_tuning), 'Folder' (default fsw_bridge/build/gains/<vehicle>).
%
%   gains.json records the controller, the law, the sizes, a SHA-256 of
%   each CSV and the git commit, so the flight software can refuse gains
%   that don't belong to its plant build (build_plant_c's plant.json has
%   the same git_commit).

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Vehicle', 'ctv');
p.addParameter('Controller', '');
p.addParameter('Folder', '');
p.parse(varargin{:});
vehicle = char(p.Results.Vehicle);
outDir = char(p.Results.Folder);
if isempty(outDir), outDir = fullfile(here, 'build', 'gains', vehicle); end
if ~isfolder(outDir), mkdir(outDir); end

if isempty(p.Results.Controller)
    V = variantConfig('VEHICLE', vehicle);
else
    V = variantConfig('VEHICLE', vehicle, 'CONTROLLER', char(p.Results.Controller));
end
A = variantAxes();
controller = A.CONTROLLER(V.VARIANT_CONTROLLER);
switch controller
    case "lqr_sched",      k1 = 'K1flat';     k2 = 'K2grid';
    case "lqr_ctv_tuning", k1 = 'K1flat_ctv'; k2 = 'K2grid_ctv';
    otherwise
        error('export_plant_gains:controller', ...
            'CONTROLLER=%s has no time-scheduled gain tables to export', controller);
end

W = hopperWorkspace([], V);
tables = {'tgrid', 'tgrid'; 'xref', 'Target_Trajectory'; 'k1flat', k1; ...
          'k2grid', k2; 'unom', 'unom'};
sizes = struct();
sha256 = struct();
for i = 1:size(tables, 1)
    T = W.(tables{i, 2});
    if isvector(T), T = T(:); end
    f = fullfile(outDir, [tables{i, 1} '.csv']);
    fid = fopen(f, 'w');
    fprintf(fid, [strjoin(repmat({'%.17g'}, 1, size(T, 2)), ',') '\n'], T.');
    fclose(fid);
    sizes.(tables{i, 1}) = size(T);
    sha256.(tables{i, 1}) = sha256File(f);
    fprintf('%-7s %-18s %d x %d\n', tables{i, 1}, tables{i, 2}, size(T, 1), size(T, 2));
end
n = sizes.tgrid(1);
assert(all(structfun(@(s) s(1), sizes) == n), 'export_plant_gains:rows', ...
    'gain tables have different row counts');

info = struct('vehicle', vehicle, 'controller', controller, ...
    'variant_controller', V.VARIANT_CONTROLLER, ...
    'gain_tables', struct('k1flat', k1, 'k2grid', k2), ...
    'law', 'u = unom - K1*(x - xref) - K2', ...
    'interpolation', 'linear in t between tgrid rows, held at both ends', ...
    'quaternion', 'scalar-first; negate x quaternion when dot(q, qref) < 0; dx is plain subtraction', ...
    'k1flat_layout', 'row-major', ...   % machine-checked by the flight software: keep exact
    'k1flat_layout_note', 'flat(13*i + j) = K1(i, j), 0-based; MATLAB reshape(K1.'', 1, 52)', ...
    'state', 'pos NED m, vel NED m/s, body rates rad/s, quaternion scalar-first (13)', ...
    'command', 'thrust N, TVC pitch rad, TVC yaw rad, RCS torque N m', ...
    'simulink_controller_rate_hz', 1000, ...
    'sizes', sizes, 'sha256', sha256, 'git_commit', gitCommit(here), ...
    'exported', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
fid = fopen(fullfile(outDir, 'gains.json'), 'w');
fprintf(fid, '%s\n', jsonencode(info, 'PrettyPrint', true));
fclose(fid);
fprintf('%s gains (%s) -> %s\n', vehicle, controller, outDir);
end


function c = gitCommit(here)
[s, r] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
if s == 0, c = strtrim(r); else, c = 'unknown'; end
[s, r] = system(sprintf('git -C "%s" status --porcelain -- ..', here));
if s == 0 && ~isempty(strtrim(r)), c = [c '-dirty']; end
end

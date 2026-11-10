function outDir = export_plant_gains(varargin)
%EXPORT_PLANT_GAINS  Write a vehicle's LQR gain schedule for the flight software.
%   export_plant_gains() writes the CTV's time-scheduled LQR tables, the
%   ones CONTROLLER=lqr_sched flies in models/hopper_6dof.slx, to
%   fsw_bridge/build/gains/ctv as CSV (one row per grid time, full double
%   precision) plus gains.json saying what they are:
%     tgrid.csv   n x 1    s, schedule times from the start of the flight
%     xref.csv    n x 13   reference state (Target_Trajectory)
%     k1flat.csv  n x 52   K1 (4 x 13, column-major) per row
%     k2grid.csv  n x 4    K2 per row
%     unom.csv    n x 4    nominal command per row
%   The harness flies u = unom - K1 * (x - xref) - K2 (LQR_Gain_Scheduler).
%
%   Options: 'Vehicle' (default 'ctv'), 'Folder' (default
%   fsw_bridge/build/gains/<vehicle>).
%
%   The workspace comes from hopperWorkspace with variantConfig('VEHICLE',
%   vehicle), the same route buildSimInput takes, so the gains match the
%   plant build_plant_c generates for that vehicle.

here = fileparts(mfilename('fullpath'));
p = inputParser;
p.addParameter('Vehicle', 'ctv');
p.addParameter('Folder', '');
p.parse(varargin{:});
vehicle = char(p.Results.Vehicle);
outDir = char(p.Results.Folder);
if isempty(outDir), outDir = fullfile(here, 'build', 'gains', vehicle); end
if ~isfolder(outDir), mkdir(outDir); end

W = hopperWorkspace([], variantConfig('VEHICLE', vehicle));
tables = {'tgrid', 'tgrid'; 'xref', 'Target_Trajectory'; 'k1flat', 'K1flat'; ...
          'k2grid', 'K2grid'; 'unom', 'unom'};
sizes = struct();
for i = 1:size(tables, 1)
    A = W.(tables{i, 2});
    if isvector(A), A = A(:); end
    fid = fopen(fullfile(outDir, [tables{i, 1} '.csv']), 'w');
    fprintf(fid, [strjoin(repmat({'%.17g'}, 1, size(A, 2)), ',') '\n'], A.');
    fclose(fid);
    sizes.(tables{i, 1}) = size(A);
    fprintf('%-7s %-18s %d x %d\n', tables{i, 1}, tables{i, 2}, size(A, 1), size(A, 2));
end

[s, r] = system(sprintf('git -C "%s" rev-parse --short HEAD', here));
if s == 0, commit = strtrim(r); else, commit = 'unknown'; end
info = struct('vehicle', vehicle, 'controller', 'lqr_sched', ...
    'law', 'u = unom - K1*(x - xref) - K2', 'state', 'pos NED, vel NED, body rates, quaternion (13)', ...
    'command', 'thrust N, TVC pitch rad, TVC yaw rad, RCS torque N m', ...
    'sizes', sizes, 'git_commit', commit, ...
    'exported', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
fid = fopen(fullfile(outDir, 'gains.json'), 'w');
fprintf(fid, '%s\n', jsonencode(info, 'PrettyPrint', true));
fclose(fid);
end

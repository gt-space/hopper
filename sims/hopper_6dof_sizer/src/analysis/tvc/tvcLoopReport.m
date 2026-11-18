function R = tvcLoopReport(file, outDir)
%TVCLOOPREPORT  Tables and figures from tvcLoopCampaign results.
%   R = tvcLoopReport()   reads results/tvc_loop.mat, prints the Monte
%   Carlo summary per config and the stress cases, and writes
%   outputs/figures/tvc/tvc_loop_*.png.

if nargin < 1, file = hopperRoot('results', 'tvc_loop.mat'); end
if nargin < 2, outDir = hopperRoot('outputs', 'figures', 'tvc'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
C = load(file);
T = C.runs;
cfg = string({C.configs.name});
mc = startsWith(T.case, ["nominal", "mc_"]);

% ---- Monte Carlo summary per config ----
rows = [];
for c = cfg
    t = T(mc & T.config == c, :);
    rows = [rows; table(c, height(t), nnz(t.success), nnz(t.meets_requirements), ...
        median(t.landing_offset, 'omitnan'), prctile(t.landing_offset, 95), max(t.max_tilt), ...
        max(t.max_cmd_deg), max(t.max_rate_dps), prctile(t.p99_rate_dps, 95), max(t.max_accel_dps2), ...
        mean(t.rms_err_deg, 'omitnan'), max(t.max_err_deg), max(t.rate_sat_pct), ...
        'VariableNames', {'config', 'runs', 'success', 'meets_req', 'offset_med_m', 'offset_p95_m', ...
        'max_tilt_deg', 'max_cmd_deg', 'max_rate_dps', 'p95_of_p99_rate_dps', 'max_accel_dps2', ...
        'mean_rms_err_deg', 'max_err_deg', 'max_rate_sat_pct'})]; %#ok<AGROW>
end
R.mc = rows;
R.stress = T(~mc, {'case', 'config', 'success', 'meets_requirements', 'landing_offset', 'max_tilt', ...
    'max_cmd_deg', 'max_ach_deg', 'max_rate_dps', 'max_accel_dps2', 'rms_err_deg', 'max_err_deg', 'rate_sat_pct', 'missed'});
R.stress = sortrows(R.stress, {'case', 'config'});
% paired difference vs the PDR model on the same scenarios
base = T(mc & T.config == cfg(1), :);
for c = cfg(2:end)
    t = T(mc & T.config == c, :);
    [~, ia, ib] = intersect(base.scenario, t.scenario);
    d = t.landing_offset(ib) - base.landing_offset(ia);
    R.paired.(matlab.lang.makeValidName(c)) = [median(d, 'omitnan'), max(abs(d))];
end
disp(R.mc); disp(R.stress);

% ---- figures ----
f = figure('Visible', 'off', 'Position', [0 0 1300 420], 'Color', 'w', 'Theme', 'light');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
short = replace(cfg, [" recommended", " sizer layout", ", 2.5 cm load", "PDR 6DOF (2nd order, 15 deg/s)"], ...
                     [" rec.", " sizer", " +2.5cm", "PDR 6DOF"]);
metrics = {'landing_offset', 'Landing offset (m)'; 'max_rate_dps', 'Peak gimbal rate (deg/s)'; 'max_err_deg', 'Peak tracking error (deg)'};
for m = 1:3
    nexttile(tl); hold on; grid on
    g = []; v = [];
    for c = 1:numel(cfg)
        t = T(mc & T.config == cfg(c), :);
        g = [g; repmat(c, height(t), 1)]; v = [v; t.(metrics{m, 1})]; %#ok<AGROW>
    end
    boxchart(g, v);
    xticks(1:numel(cfg)); xticklabels(short); xtickangle(30);
    ylabel(metrics{m, 2});
end
title(tl, sprintf('Model-in-the-loop Monte Carlo (%d dispersed + nominal per config)', C.n));
exportgraphics(f, fullfile(outDir, 'tvc_loop_mc.png'), 'Resolution', 200); close(f)
end

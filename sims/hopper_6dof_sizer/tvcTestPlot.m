%% Plot TVC clamp campaign results
if ~exist('C','var')
    S = load(fullfile('results','tvc_clamp_mc.mat'));
    if isfield(S,'C'), C = S.C; else, C = S; end
end
R = C.runs;
cases = unique(R.case, 'stable');
cfgs  = unique(R.config, 'stable');

getM = @(col) gridOf(R, col, cases, cfgs);
caseLbl = strrep(cases, '_', ' ');

figure('Name','TVC clamp campaign','Position',[100 100 1300 750])
tiledlayout(2,3,'TileSpacing','compact')

nexttile; bar(getM('max_tilt')); set(gca,'YScale','log');
title('Max tilt (deg)'); grid on; xticklabels(caseLbl)

nexttile; bar(getM('rms_tilt_deg'));
title('RMS tilt (deg)'); grid on; xticklabels(caseLbl)

nexttile; bar(getM('landing_offset'));
title('Landing offset (m)'); grid on; xticklabels(caseLbl)

nexttile; bar(getM('at_rate_lim_pct'));
title('% of flight at rate clamp'); grid on; xticklabels(caseLbl)

nexttile; bar(getM('accel_over_lim')); yline(1,'--','clamp');
title('Peak accel / accel clamp'); grid on; xticklabels(caseLbl)

nexttile; bar(getM('osc_amp_deg')); yline(0.1,'--','osc threshold');
set(gca,'YScale','log'); title('Gimbal oscillation amplitude (deg)'); grid on; xticklabels(caseLbl)

lg = legend(cfgs, 'Orientation','horizontal');  lg.Layout.Tile = 'south';
set(findall(gcf,'Type','axes'), 'TickLabelInterpreter','none', 'XTickLabelRotation', 20)

function M = gridOf(R, col, cases, cfgs)
M = nan(numel(cases), numel(cfgs));
for i = 1:numel(cases)
    for j = 1:numel(cfgs)
        v = R.(col)(R.case == cases(i) & R.config == cfgs(j));
        if ~isempty(v), M(i,j) = v(1); end
    end
end
end
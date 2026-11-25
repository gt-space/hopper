function R = tvcDfSweep(opts)
%TVCDFSWSWEEP  Describing-function sweep: gain + phase vs. amplitude × frequency.
%   R = tvcDfSweep()
%
%   Runs the full E+L electromechanical model (tvcAxisSim) at every
%   (amplitude, frequency) grid point for each actuator configuration.
%   For each point, the steady-state output is projected onto sin/cos at
%   the drive frequency to extract gain and phase (same method as
%   tvcTrade.sineResp). The drive acceleration ramp is set to TVC team's
%   values via tvcAxisSim AccelLimit.
%
%   Key outputs (all in deg / deg/s / deg/s^2):
%     R.pts   – table of all (config, A_deg, f_Hz, gain_dB, phase_deg)
%     R.fig_phase_vs_A  – figure: phase at 2 Hz vs. amplitude, all configs
%     R.fig_bode        – figure: gain & phase vs. freq, A = 0.5 / 2°
%
%   Configs
%     AG 9113, accel limit 235 deg/s^2  (TVC team number, ±8° layout)
%     AG 9313, accel limit 640 deg/s^2  (TVC team number, ±8° layout)
%     AG 9113, unlimited                (physical motor limits only)
%     GNC floor, accel limit 200 deg/s^2 (minimum below which 2 Hz tracking fails)

arguments
    opts.OutDir  (1,:) char = fullfile(hopperRootSafe(), 'outputs', 'figures', 'tvc')
    opts.SaveFile (1,:) char = fullfile(hopperRootSafe(), 'results', 'tvc_df_sweep.mat')
    opts.Amps_deg (1,:) double = [0.25 0.5 1 1.5 2 3]
    opts.Freqs_Hz (1,:) double = [0.5 1 2 3 5]
end

if ~exist(opts.OutDir, 'dir'), mkdir(opts.OutDir); end

REQ = tvcRequirements();

% ---- Configurations ----------------------------------------------------------------
%  name, actuator, layout, RateLimit (deg/s), AccelLimit (deg/s^2)
cfgs = {
    'AG 9113 | 20.5°/s, 235°/s²',   'AG9113', 'recommended', 20.5, 235
    'AG 9313 | 13.5°/s, 640°/s²',   'AG9313', 'recommended', 13.5, 640
    'AG 9113 | unlimited',            'AG9113', 'recommended',  Inf, Inf
    'GNC floor | 25°/s, 200°/s²',    'AG9113', 'recommended', 25.0, 200
};
nC = size(cfgs, 1);
A_deg = opts.Amps_deg;
f_Hz  = opts.Freqs_Hz;
nA = numel(A_deg);  nF = numel(f_Hz);

% ---- Sweep -------------------------------------------------------------------------
gain_dB   = nan(nC, nA, nF);
phase_deg = nan(nC, nA, nF);

for c = 1:nC
    P = tvcActuatorParams(cfgs{c,2}, Layout=cfgs{c,3});
    rlim = deg2rad(cfgs{c,4});
    alim = deg2rad(cfgs{c,5});
    fprintf('Config %d/%d: %s  (r_eff=%.1f mm, rate_noload=%.1f °/s)\n', ...
        c, nC, cfgs{c,1}, P.r_eff*1e3, rad2deg(P.rate_noload));
    for fi = 1:nF
        f = f_Hz(fi);
        for ai = 1:nA
            A = A_deg(ai);
            r = sinePoint(P, A, f, rlim, alim);
            gain_dB(c, ai, fi)   = r.gain_dB;
            phase_deg(c, ai, fi) = r.phase_deg;
        end
    end
end

% ---- Pack results ------------------------------------------------------------------
rows = {};
for c = 1:nC
    for fi = 1:nF
        for ai = 1:nA
            rows{end+1} = {cfgs{c,1}, A_deg(ai), f_Hz(fi), gain_dB(c,ai,fi), phase_deg(c,ai,fi)}; %#ok<AGROW>
        end
    end
end
R.pts = cell2table(vertcat(rows{:}), 'VariableNames', {'config','A_deg','f_Hz','gain_dB','phase_deg'});
R.configs  = cfgs;
R.A_deg    = A_deg;
R.f_Hz     = f_Hz;
R.gain_dB  = gain_dB;
R.phase_deg = phase_deg;
R.REQ      = REQ;

% ---- Print key table ---------------------------------------------------------------
f2 = find(f_Hz == 2);
fprintf('\n==== DF at 2 Hz ====\n');
fprintf('%-28s  %5s  %6s  %7s\n', 'Config', 'A(°)', 'Gain dB', 'Phase °');
for c = 1:nC
    for ai = 1:nA
        fprintf('%-28s  %5.2f  %6.2f  %7.1f\n', cfgs{c,1}, A_deg(ai), gain_dB(c,ai,f2), phase_deg(c,ai,f2));
    end
end

% ---- Find limiting amplitude (phase < -150 at 2 Hz) --------------------------------
fprintf('\n==== Phase < -150° at 2 Hz ====\n');
for c = 1:nC
    ph = phase_deg(c, :, f2);
    idx = find(ph < -150, 1);
    if isempty(idx)
        fprintf('  %s: phase stays above -150° across all amplitudes (min %.0f° at %.2f°)\n', cfgs{c,1}, min(ph), A_deg(end));
    else
        fprintf('  %s: phase drops below -150° at A >= %.2f° (%.0f°)\n', cfgs{c,1}, A_deg(idx), ph(idx));
    end
end

% ---- Figures -----------------------------------------------------------------------
clrs = [0.13 0.27 0.53; 0.85 0.55 0.0; 0.40 0.70 0.40; 0.65 0.15 0.15];
ls   = {'-', '-', '--', ':'};

% Figure 1: phase vs amplitude at 2 Hz
fig1 = figure('Visible','off','Position',[50 50 700 420],'Color','w','Theme','light');
ax = axes(fig1);  hold(ax,'on');  grid(ax,'on');
yline(ax, -150, 'k:', 'LineWidth', 1.2, 'Label', '-150° limit');
yline(ax, REQ.phase2Hz_deg, 'k--', 'LineWidth', 1, 'Label', sprintf('req %d°', REQ.phase2Hz_deg));
for c = 1:nC
    plot(ax, A_deg, squeeze(phase_deg(c,:,f2)), [ls{c},'o'], 'Color', clrs(c,:), 'LineWidth', 1.6, ...
        'DisplayName', cfgs{c,1});
end
xlabel(ax, 'Command amplitude (deg)');
ylabel(ax, 'Phase at 2 Hz (deg)');
title(ax, 'TVC DF sweep — phase vs. amplitude at 2 Hz');
legend(ax, 'Location', 'southwest');
exportgraphics(fig1, fullfile(opts.OutDir,'tvc_df_phase_vs_amp.png'), 'Resolution', 200);
R.fig_phase_vs_A = fig1;

% Figure 2: Bode-like, A = 0.5° and 2°
a_bode = [0.5 2];  lsB = {'-', '--'};
fig2 = figure('Visible','off','Position',[50 50 800 600],'Color','w','Theme','light');
tl = tiledlayout(fig2, 2, 1, 'TileSpacing','compact');

ax1 = nexttile(tl);  hold(ax1,'on');  grid(ax1,'on');  set(ax1,'XScale','log');
ax2 = nexttile(tl);  hold(ax2,'on');  grid(ax2,'on');  set(ax2,'XScale','log');
yline(ax1, -1, 'k:', 'Label', '-1 dB');
yline(ax2, -150, 'k:', 'Label', '-150°');
yline(ax2, REQ.phase2Hz_deg, 'k--', 'Label', sprintf('req %d°', REQ.phase2Hz_deg));

for c = 1:nC
    for bi = 1:numel(a_bode)
        ai = find(A_deg == a_bode(bi));
        if isempty(ai), continue; end
        lbl = sprintf('%s, A=%.1f°', cfgs{c,1}, a_bode(bi));
        col = clrs(c,:) * (0.5 + 0.5*bi/numel(a_bode));
        semilogx(ax1, f_Hz, squeeze(gain_dB(c,ai,:)), [lsB{bi},'o'], 'Color', col, 'LineWidth', 1.4, 'DisplayName', lbl);
        semilogx(ax2, f_Hz, squeeze(phase_deg(c,ai,:)), [lsB{bi},'o'], 'Color', col, 'LineWidth', 1.4, 'HandleVisibility','off');
    end
end
xline(ax1, REQ.wc_rad_s/(2*pi), ':', 'pitch Xover', 'HandleVisibility','off');
xline(ax2, REQ.wc_rad_s/(2*pi), ':', 'pitch Xover');
ylabel(ax1, 'Gain (dB)');  ylabel(ax2, 'Phase (deg)');
xlabel(ax2, 'Frequency (Hz)');  title(ax1, 'TVC DF sweep — frequency response');
legend(ax1, 'Location','southwest','FontSize',7);
exportgraphics(fig2, fullfile(opts.OutDir,'tvc_df_bode.png'), 'Resolution', 200);
R.fig_bode = fig2;

save(opts.SaveFile, '-struct', 'R', '-v7.3');
fprintf('\nSaved %s\n', opts.SaveFile);
fprintf('Figures: %s\n', opts.OutDir);
end

%% ---- helpers -----------------------------------------------------------------------
function r = sinePoint(P, A_deg, f_Hz, rate_lim_rads, accel_lim_rads2)
% Run tvcAxisSim with a sine command and extract gain + phase.
A = deg2rad(A_deg);  f = f_Hz;  w = 2*pi*f;
nc = max(5, ceil(f * 1.5));     % steady-state cycles to keep
t_end = 0.3 + (nc + 3) / f;    % settle + nc full periods
S = tvcAxisSim(P, [], t_end, @(t) A * sin(w * t), accel_lim=accel_lim_rads2, rate_lim=rate_lim_rads);
m = S.t > t_end - nc/f;
X = [sin(w * S.t(m)), cos(w * S.t(m))];
c = X \ S.delta(m);
r.gain_dB   = 20 * log10(hypot(c(1), c(2)) / A);
r.phase_deg = rad2deg(atan2(c(2), c(1)));
end

function r = hopperRootSafe()
try
    r = hopperRoot();
catch
    r = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
end

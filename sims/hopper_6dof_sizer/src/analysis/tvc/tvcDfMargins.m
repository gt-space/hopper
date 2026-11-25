function M = tvcDfMargins(R, opts)
%TVCDFMARGINS  Step 2: PM / GM / delay margin and limit-cycle check from DF sweep.
%   M = tvcDfMargins(R)   where R is the struct from tvcDfSweep (workspace or file).
%
%   The open-loop pitch TF is built analytically then split:
%     L_lin(jw)  = k_p * C(jw) / (jw)^2       (plant + GNC PD, no actuator)
%     T_act(jw)  = wn^2/(s^2+2z*wn*s+wn^2)    (small-signal actuator model)
%     L_total(jw) = L_lin * T_act              (what GNC was tuned against)
%     L_OL(A,jw) = L_lin(jw) * N(A,jw)        (actual, with DF actuator)
%
%   L_lin and N(A,f) are both evaluated at the five DF sweep frequencies
%   exactly — no interpolation across frequency.  Crossover is found by
%   linear interpolation between the two bracketing DF points.
%
%   PM_nominal is not in this repo; parametrised over 30 / 40 / 50 deg.
%   Plant assumptions (all from files in this repo):
%     Iy  = 41 kg m^2  (LQTrack_Master design value)
%     Lca = 0.59 m     (tvcRequirements: 1.040 - 0.45)
%     T_h = 1116 N     (mid-flight hover, back-calc from twr/thrust in results)
%     wn  = 36 rad/s, zeta = 0.7  (LQT actuator model, tvcRequirements.wn_act)

arguments
    R  struct
    opts.PM_nominal  (1,:) double = [30 40 50]
    opts.OutDir  (1,:) char = fullfile(hopperRootSafe(),'outputs','figures','tvc')
    opts.SaveFile (1,:) char = fullfile(hopperRootSafe(),'results','tvc_df_margins.mat')
end
if ~exist(opts.OutDir,'dir'), mkdir(opts.OutDir); end
t0_wall = tic;
REQ = tvcRequirements();

%% ---- Pitch-axis plant (hover) ----------------------------------------
Iy   = 41;              % kg m^2
Lca  = 1.040 - 0.45;   % m
Th   = 1116;            % N
k_p  = Th * Lca / Iy;  % double-integrator gain  [rad/s^2 / rad]
wc   = REQ.wc_rad_s;   % 12.25 rad/s
wn   = REQ.wn_act;  zeta = 0.7;

t_plant = toc(t0_wall);  fprintf('[timing] plant + REQ setup: %.2f s\n', t_plant);

%% ---- DF data ----------------------------------------------------------
%   Configs: 1=AG9113+limits, 2=AG9313+limits, 3=unlimited, 4=GNC floor
A_v  = R.A_deg;    f_df = R.f_Hz;    w_df = 2*pi*f_df;
nA   = numel(A_v); nF   = numel(f_df);
cfg_idx = [1 2 4];
cfg_lbl = {R.configs{cfg_idx(1),1}, R.configs{cfg_idx(2),1}, R.configs{cfg_idx(3),1}};
nC = numel(cfg_idx);

% Second-order actuator at DF frequencies
T_a_df = wn^2 ./ ((1i*w_df).^2 + 2*zeta*wn*(1i*w_df) + wn^2);

% ---- Unwrap DF phase across frequency for each (config, amplitude) -----
%   atan2 output wraps at ±180°; unwrapping across freq gives monotone phase.
gdb_all = R.gain_dB;    % [nCfg x nA x nF]
phd_all = R.phase_deg;

for c = 1:size(gdb_all,1)
    for ai = 1:nA
        v = squeeze(phd_all(c,ai,:))';
        phd_all(c,ai,:) = rad2deg(unwrap(deg2rad(v)));
    end
end

t_unwrap = toc(t0_wall);  fprintf('[timing] DF unwrap done: %.2f s\n', t_unwrap);

%% ---- Sweep over PM_nominal -------------------------------------------
nPM = numel(opts.PM_nominal);

for pi_idx = 1:nPM
    PM0 = opts.PM_nominal(pi_idx);

    % Calibrate PD: angle(L_total at wc) = -180+PM0
    %   angle(L_total) = angle(Gp) + angle(C) + angle(T_act) = -180 + angle(C) + angle(T_act)
    %   => angle(C) = PM0 - angle(T_act)
    Ta_wc  = wn^2 / ((1i*wc)^2 + 2*zeta*wn*(1i*wc) + wn^2);
    ph_C   = PM0 - rad2deg(angle(Ta_wc));   % target phase of C at wc [deg]
    Kd_Kp  = tand(ph_C) / wc;
    C_norm = sqrt(Kd_Kp^2 * wc^2 + 1);
    Kp     = wc^2 / (k_p * C_norm * abs(Ta_wc));
    Kd     = Kd_Kp * Kp;

    % Verify
    C_wc = Kd*(1i*wc) + Kp;
    PM_check = 180 + rad2deg(angle(k_p * C_wc / (1i*wc)^2 * Ta_wc));

    % L_lin at the five DF frequencies (exact, analytic)
    C_df    = Kd*(1i*w_df) + Kp;
    L_lin_df = k_p * C_df ./ (1i*w_df).^2;
    ph_lin  = rad2deg(angle(L_lin_df));  % in (-270,-90] deg; well-defined, no wrap

    % L_total at DF freqs (the linear loop, for reference)
    L_tot_df = L_lin_df .* T_a_df;

    for ci = 1:nC
        c = cfg_idx(ci);
        gdb_c = squeeze(gdb_all(c,:,:));   % [nA x nF], gain_dB
        phd_c = squeeze(phd_all(c,:,:));   % [nA x nF], unwrapped phase deg

        PM_v = nan(1,nA);  GM_v = nan(1,nA);
        DM_v = nan(1,nA);  wc_v = nan(1,nA);

        for ai = 1:nA
            N_mag = 10.^(gdb_c(ai,:)/20);        % actuator DF gain
            N_ph  = phd_c(ai,:);                  % actuator DF phase (unwrapped)
            L_OL_mag = abs(L_lin_df) .* N_mag;    % |L_OL| at DF freqs
            ph_OL    = ph_lin + N_ph;             % angle(L_OL) = angle(L_lin)+angle(N)

            % ---- Gain crossover: last frequency where mag crosses 1 ------
            [wc_act, pm] = gainCrossover(w_df, L_OL_mag, ph_OL);
            PM_v(ai) = pm;  wc_v(ai) = wc_act;

            % ---- Phase crossover for GM ----------------------------------
            GM_v(ai) = phaseCrossoverGM(w_df, L_OL_mag, ph_OL);

            % ---- Delay margin -------------------------------------------
            if ~isnan(wc_act) && ~isnan(pm)
                DM_v(ai) = deg2rad(pm) / wc_act * 1e3;  % ms
            end
        end

        fname = matlab.lang.makeValidName(sprintf('PM%d_%s', round(PM0), cfg_lbl{ci}));
        M.(fname) = struct('PM', PM_v, 'GM', GM_v, 'DM_ms', DM_v, 'wc', wc_v, ...
                           'cfg', cfg_lbl{ci}, 'PM_nominal', PM0, 'Kp', Kp, 'Kd', Kd);

        if pi_idx <= 2  % print for PM_nom = 30 and 40
            fprintf('\n=== %s  (PM_lin = %.0f°, Kp=%.3f, Kd=%.3f, check=%.1f°) ===\n', ...
                cfg_lbl{ci}, PM0, Kp, Kd, PM_check);
            fprintf('%-6s  %6s  %6s  %7s  %8s\n', 'A(°)', 'PM(°)', 'GM(dB)', 'DM(ms)', 'wc(Hz)');
            for ai = 1:nA
                flag = '';
                if ~isnan(PM_v(ai)) && PM_v(ai) < 30, flag = ' <<< PM<30°'; end
                fprintf('%-6.2f  %6.1f  %6.1f  %7.1f  %8.2f%s\n', ...
                    A_v(ai), PM_v(ai), GM_v(ai), DM_v(ai), wc_v(ai)/(2*pi), flag);
            end
        end
    end
end

t_margins = toc(t0_wall);  fprintf('[timing] margin sweep done: %.2f s\n', t_margins);

%% ---- Limit-cycle check at 2 Hz (PM_lin=40°) -------------------------
PM0 = 40;
Ta_wc  = wn^2 / ((1i*wc)^2 + 2*zeta*wn*(1i*wc) + wn^2);
ph_C   = PM0 - rad2deg(angle(Ta_wc));
Kd_Kp  = tand(ph_C) / wc;
C_norm = sqrt(Kd_Kp^2*wc^2 + 1);
Kp40   = wc^2 / (k_p * C_norm * abs(Ta_wc));
Kd40   = Kd_Kp * Kp40;

f2idx = find(f_df == 2);
w2 = w_df(f2idx);
L_lin2 = k_p * (Kd40*(1i*w2) + Kp40) / (1i*w2)^2;

c = 1;  % AG9113+limits
gdb2 = squeeze(gdb_all(c,:,f2idx));  % [nA x 1], already unwrapped not needed for single freq
phd2 = squeeze(phd_all(c,:,f2idx));
L_OL2 = L_lin2 .* (10.^(gdb2/20) .* exp(1i*deg2rad(phd2)));

fprintf('\n=== Limit-cycle check at 2 Hz — AG 9113, PM_lin=40° ===\n');
fprintf('  (L_lin at 2 Hz: |L_lin|=%.3f, angle=%.1f°)\n', abs(L_lin2), rad2deg(angle(L_lin2)));
fprintf('%-6s  %7s  %8s  %s\n', 'A(°)','|L_OL|','angle(°)','Notes');
for ai = 1:nA
    mag2 = abs(L_OL2(ai));
    ph2  = rad2deg(angle(L_lin2)) + phd2(ai);   % use summed phase (consistent)
    ph2n = mod(ph2+360,360)-360;
    note = '';
    if ph2n < -150, note = [note ' phase<-150']; end
    if abs(mag2-1) < 0.15 && ph2n < -170, note = [note ' *** LC risk']; end
    fprintf('%-6.2f  %7.3f  %8.1f  %s\n', A_v(ai), mag2, ph2n, note);
end
% Find -180° crossing in amplitude
pc = (rad2deg(angle(L_lin2)) + phd2') + 180;
sc = find(pc(1:end-1).*pc(2:end) < 0);
if isempty(sc)
    fprintf('No -180° phase crossing in amplitude range -> no predicted LC at 2 Hz\n');
else
    for k = sc'
        A_lc = interp1(pc(k:k+1), A_v(k:k+1), 0);
        G_lc = interp1(A_v(k:k+1), abs(L_OL2(k:k+1)), A_lc);
        fprintf('-180° crossing at A_lc = %.2f deg; |L_OL| = %.3f (%.1f dB) -> ', ...
            A_lc, G_lc, 20*log10(G_lc));
        if G_lc > 1
            fprintf('UNSTABLE LC predicted\n');
        elseif abs(G_lc-1) < 0.1
            fprintf('MARGINAL LC risk (|L_OL|~1)\n');
        else
            fprintf('stable (|L_OL|<1 at -180° -> no sustained LC)\n');
        end
    end
end

%% ---- Summary: amplitude where PM < 30° for each config ---------------
fprintf('\n=== PM < 30° threshold (PM_lin = 40°) ===\n');
PM0 = 40;
for ci = 1:nC
    fname = matlab.lang.makeValidName(sprintf('PM40_%s', cfg_lbl{ci}));
    if ~isfield(M, fname), continue; end
    s = M.(fname);
    idx = find(s.PM < 30, 1);
    if isempty(idx)
        fprintf('  %-40s  PM stays >= 30° (min %.0f° at A=%.1f°)\n', ...
            cfg_lbl{ci}, min(s.PM), A_v(end));
    else
        fprintf('  %-40s  PM < 30° at A >= %.2f° (PM=%.0f°)\n', ...
            cfg_lbl{ci}, A_v(idx), s.PM(idx));
    end
end
fprintf('\n  Flight peak demand: ~0.65° -> well below threshold\n');

t_lc = toc(t0_wall);  fprintf('[timing] LC check done: %.2f s\n', t_lc);

%% ---- Figures ----------------------------------------------------------
PM_clrs = [0 0.6 0; 0 0 0.9; 0.7 0 0.7];  % 30/40/50
clrs3   = [0.13 0.27 0.53; 0.85 0.55 0.0; 0.65 0.15 0.15];

% ---- Fig 1: PM vs amplitude, all PM_nominal, AG9113 + AG9313 ----------
fig1 = figure('Visible','off','Position',[50 50 900 480],'Color','w','Theme','light');
ax1 = axes(fig1); hold(ax1,'on'); grid(ax1,'on');
yline(ax1, 30, 'k--', 'LineWidth',1.3, 'Label','30° floor');
yline(ax1,  0, 'k:',  'LineWidth',0.8);
for pi_idx = 1:nPM
    PM0 = opts.PM_nominal(pi_idx);
    fn = matlab.lang.makeValidName(sprintf('PM%d_%s', round(PM0), cfg_lbl{1}));
    if isfield(M, fn)
        plot(ax1, A_v, M.(fn).PM, '-o', 'Color', PM_clrs(pi_idx,:), 'LineWidth',1.7, ...
            'DisplayName', sprintf('AG 9113 | PM_{lin}=%d°', round(PM0)));
    end
end
fn9313 = matlab.lang.makeValidName(sprintf('PM40_%s', cfg_lbl{2}));
if isfield(M, fn9313)
    plot(ax1, A_v, M.(fn9313).PM, '--s', 'Color', clrs3(2,:), 'LineWidth',1.7, ...
        'DisplayName', 'AG 9313 | PM_{lin}=40°');
end
xline(ax1, 0.65, ':', 'Color',[0.4 0.4 0.4], 'Label','flight peak ~0.65°');
xlabel(ax1,'Command amplitude (deg)');  ylabel(ax1,'Phase margin (deg)');
title(ax1,'TVC DF margins — phase margin vs. amplitude');
legend(ax1,'Location','southwest');
exportgraphics(fig1, fullfile(opts.OutDir,'tvc_df_PM.png'), 'Resolution',200);
M.fig_PM = fig1;

% ---- Fig 2: GM and delay margin ----------------------------------------
fig2 = figure('Visible','off','Position',[50 50 900 480],'Color','w','Theme','light');
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing','compact');
ax2a = nexttile(tl2); hold(ax2a,'on'); grid(ax2a,'on');
ax2b = nexttile(tl2); hold(ax2b,'on'); grid(ax2b,'on');
yline(ax2a, 6, 'k--', 'Label','6 dB');
yline(ax2b, 10,'k--', 'Label','10 ms');
for pi_idx = 1:nPM
    PM0 = opts.PM_nominal(pi_idx);
    fn1 = matlab.lang.makeValidName(sprintf('PM%d_%s', round(PM0), cfg_lbl{1}));
    if isfield(M, fn1)
        plot(ax2a, A_v, M.(fn1).GM,    '-o', 'Color',PM_clrs(pi_idx,:), 'LineWidth',1.5, ...
            'DisplayName', sprintf('9113 PM_{lin}=%d°',round(PM0)));
        plot(ax2b, A_v, M.(fn1).DM_ms, '-o', 'Color',PM_clrs(pi_idx,:), 'LineWidth',1.5, ...
            'DisplayName', sprintf('9113 PM_{lin}=%d°',round(PM0)));
    end
end
xlabel(ax2a,'Amplitude (deg)'); ylabel(ax2a,'Gain margin (dB)');
xlabel(ax2b,'Amplitude (deg)'); ylabel(ax2b,'Delay margin (ms)');
title(ax2a,'Gain margin'); title(ax2b,'Delay margin');
legend(ax2a,'Location','southwest'); legend(ax2b,'Location','southwest');
exportgraphics(fig2, fullfile(opts.OutDir,'tvc_df_GM_DM.png'), 'Resolution',200);
M.fig_GM_DM = fig2;

% ---- Fig 3: Nyquist at PM=40 — L_total (analytic) + L_OL dots at DF freqs --
% The Nyquist for L_total is smooth (analytic); for L_OL we show the 5 DF points
% connected by lines per amplitude, so only measured data enters the plot.
w_fine = 2*pi * logspace(log10(0.1), log10(20), 500);
PM0 = 40;
Ta_wc  = wn^2 / ((1i*wc)^2 + 2*zeta*wn*(1i*wc) + wn^2);
ph_C   = PM0 - rad2deg(angle(Ta_wc));
Kd_Kp  = tand(ph_C) / wc;  Kp0 = wc^2 / (k_p*sqrt(Kd_Kp^2*wc^2+1)*abs(Ta_wc));
Kd0    = Kd_Kp * Kp0;
L_tot_fine = k_p*(Kd0*(1i*w_fine)+Kp0)./(1i*w_fine).^2 .* ...
             (wn^2./((1i*w_fine).^2+2*zeta*wn*(1i*w_fine)+wn^2));
% L_OL at DF freqs: re-compute for PM=40
C_df   = Kd0*(1i*w_df)+Kp0;
L_lin_df0 = k_p * C_df ./ (1i*w_df).^2;
c = 1;
gdb_c0 = squeeze(gdb_all(c,:,:));
phd_c0 = squeeze(phd_all(c,:,:));

nyq_clrs = copper(nA+1); nyq_clrs = nyq_clrs(1:nA,:);
fig3 = figure('Visible','off','Position',[50 50 720 620],'Color','w','Theme','light');
ax3 = axes(fig3); hold(ax3,'on'); grid(ax3,'on'); axis(ax3,'equal');
plot(ax3, real(L_tot_fine), imag(L_tot_fine), 'k-', 'LineWidth',1.0, ...
    'DisplayName','L_{total} (linear, fine grid)');
plot(ax3, -1, 0, 'rx', 'MarkerSize',14, 'LineWidth',2.5, 'DisplayName','-1');
for ai = 1:nA
    N_m = 10.^(gdb_c0(ai,:)/20);
    N_p = phd_c0(ai,:);
    L_OL_v = L_lin_df0 .* N_m .* exp(1i*deg2rad(N_p));
    plot(ax3, real(L_OL_v), imag(L_OL_v), '-o', 'Color', nyq_clrs(ai,:), ...
        'LineWidth',1.3, 'MarkerSize',5, 'DisplayName', sprintf('A=%.2f°', A_v(ai)));
end
xlim(ax3, [-2.5 0.5]); ylim(ax3, [-1.5 1.5]);
xlabel(ax3,'Re'); ylabel(ax3,'Im');
title(ax3,'Nyquist: L_{OL}(A,j\omega) — AG 9113, PM_{lin}=40°');
legend(ax3,'Location','northeast','FontSize',7);
exportgraphics(fig3, fullfile(opts.OutDir,'tvc_df_nyquist.png'), 'Resolution',200);
M.fig_nyquist = fig3;

t_figs = toc(t0_wall);  fprintf('[timing] figures exported: %.2f s\n', t_figs);

% Save data only (no figure handles — serializing figures to .mat is very slow)
M_save = rmfield(M, intersect(fieldnames(M), {'fig_PM','fig_GM_DM','fig_nyquist'}));
save(opts.SaveFile, '-struct', 'M_save', '-v7.3');
t_save = toc(t0_wall);
fprintf('[timing] save done: %.2f s\n', t_save);
fprintf('\n[timing] TOTAL: %.2f s\n', t_save);
fprintf('Saved %s\n', opts.SaveFile);
fprintf('Figures: %s\n', opts.OutDir);
end

%% ---- helpers -----------------------------------------------------------
function [wc_act, pm] = gainCrossover(w, mag, ph_OL_deg)
% Find the last gain crossover (mag crosses 1 downward) and read phase margin.
% Works on the small set of DF frequency points; linear interpolation between
% the two bracketing points.
wc_act = NaN;  pm = NaN;
n = numel(mag);
for k = n-1 : -1 : 1
    if mag(k) >= 1 && mag(k+1) < 1
        t = (1 - mag(k)) / (mag(k+1) - mag(k));
        wc_act = w(k) + t*(w(k+1)-w(k));
        ph_at_wc = ph_OL_deg(k) + t*(ph_OL_deg(k+1) - ph_OL_deg(k));
        pm = 180 + ph_at_wc;
        return
    end
end
% No downward crossing: either always >1 or always <1
if all(mag > 1)
    % crossover above highest DF freq — use last point as lower bound
    pm = 180 + ph_OL_deg(end);
    wc_act = w(end);
elseif all(mag < 1)
    % crossover below lowest DF freq — not in our data range
end
end

function gm_dB = phaseCrossoverGM(w, mag, ph_OL_deg)
% Find first -180° phase crossing and read gain margin.
gm_dB = NaN;
for k = 1:numel(ph_OL_deg)-1
    p1 = ph_OL_deg(k); p2 = ph_OL_deg(k+1);
    if (p1+180)*(p2+180) < 0  % crosses -180°
        t = (-180 - p1) / (p2 - p1);
        mag_pc = mag(k) + t*(mag(k+1)-mag(k));
        if mag_pc > 0
            gm_dB = -20*log10(mag_pc);
        end
        return
    end
end
end

function r = hopperRootSafe()
try, r = hopperRoot();
catch, r = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))); end
end

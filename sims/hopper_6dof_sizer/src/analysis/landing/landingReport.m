function Sm = landingReport(T, opts)
%LANDINGREPORT  Summarize a landingSweep: margins, stroke and g per core.
%   Sm = landingReport()      loads results/landing_sweep.mat
%   Sm = landingReport(T)
%
%   One row per crush-core configuration, worst case over surfaces:
%     nominal_crush_mm   largest core/pad stroke in the nominal envelope
%                        (> 0.5 mm means the core is used up in a normal
%                        landing and must be replaced)
%     off_*              off-nominal (10 ft/s) worst values
%     MS_main_buckle     Pcr / (FS_buckle * N_main) - 1
%     MS_aux_buckle      same for the aux struts in compression
%     MS_aux_yield       Sy A / (FS_yield * N_aux tension) - 1
%     core_length_mm     stroke used / 0.75 (usable fraction)
%   Writes outputs/landing/landing_sweep.csv, landing_summary.csv and
%   figures.

arguments
    T table = table()
    opts.FS_buckle (1,1) double = 2.0     % leg hand calcs
    opts.FS_yield (1,1) double = 1.5
    opts.OutDir (1,:) char = fullfile(hopperRoot(), 'outputs', 'landing')
end
if isempty(T)
    L = load(fullfile(hopperRoot(), 'results', 'landing_sweep.mat'), 'T');
    T = L.T;
end
if ~exist(opts.OutDir, 'dir'), mkdir(opts.OutDir); end
G = landingGear();
Pcr_main = G.Pcr(1, 1) / 1e3;  Pcr_aux = G.Pcr(2, 1) / 1e3;      % kN
Ty_aux = G.sy * G.A(2) / 1e3;  Cy_main = G.sy * G.A(1) / 1e3;

cfg = unique(T(:, {'placement', 'plateau_kN'}), 'rows', 'stable');
Sm = cfg;
for r = 1:height(cfg)
    m = T.placement == cfg.placement(r) & T.plateau_kN == cfg.plateau_kN(r);
    nom = T(m & T.set == "nominal", :);
    off = T(m & T.set == "offnominal", :);
    Sm.nominal_crush_mm(r) = max([nom.stroke_core_mm; nom.stroke_pad_mm]);
    Sm.nom_main_kN(r) = max(nom.main_comp_kN);
    Sm.nom_aux_tens_kN(r) = max(nom.aux_tens_kN);
    Sm.nom_acc_g(r) = max(nom.acc_lp_g);
    Sm.off_main_kN(r) = max(off.main_comp_kN);
    Sm.off_aux_tens_kN(r) = max(off.aux_tens_kN);
    Sm.off_aux_comp_kN(r) = max(off.aux_comp_kN);
    Sm.off_foot_Fn_kN(r) = max(off.foot_Fn_lp_kN);
    Sm.off_foot_Ft_kN(r) = max(off.foot_Ft_kN);
    Sm.off_acc_g(r) = max(off.acc_lp_g);
    used = max([off.stroke_core_mm; off.stroke_pad_mm]);
    Sm.off_stroke_mm(r) = used;
    Sm.core_length_mm(r) = used / 0.75;
    Sm.off_tilt_deg(r) = max(off.tilt_max_deg);
    Sm.tipped(r) = any(T.tipped(m));
    worst_main = max([nom.main_comp_kN; off.main_comp_kN]);
    worst_auxc = max([nom.aux_comp_kN; off.aux_comp_kN]);
    worst_auxt = max([nom.aux_tens_kN; off.aux_tens_kN]);
    Sm.MS_main_buckle(r) = Pcr_main / (opts.FS_buckle * worst_main) - 1;
    Sm.MS_aux_buckle(r) = Pcr_aux / (opts.FS_buckle * max(worst_auxc, eps)) - 1;
    Sm.MS_aux_yield(r) = Ty_aux / (opts.FS_yield * worst_auxt) - 1;
    Sm.MS_main_yield(r) = Cy_main / (opts.FS_yield * worst_main) - 1;
end
Sm.all_margins_ok = Sm.MS_main_buckle >= 0 & Sm.MS_aux_buckle >= 0 & Sm.MS_aux_yield >= 0 & ~Sm.tipped;
Sm.reusable = Sm.nominal_crush_mm < 0.5;

writetable(T, fullfile(opts.OutDir, 'landing_sweep.csv'));
writetable(Sm, fullfile(opts.OutDir, 'landing_summary.csv'));

% ---- Figures ----
surf = unique(T.surface, 'stable');
place = ["main", "pad"];
f = figure('Visible', 'off', 'Position', [0 0 1500 900], 'Color', 'w');
try, theme(f, 'light'); catch, end      % readable when the desktop is in dark mode
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact');
title(tl, sprintf('Crush-core sweep, %s legs: off-nominal worst case per surface (dashed: nominal)', G.joints));
for p = 1:2
    for q = 1:3
        ax = nexttile(tl, (p-1)*3 + q);  hold(ax, 'on');  grid(ax, 'on');
        for s = 1:numel(surf)
            m = T.placement == place(p) & T.surface == surf(s);
            F = unique(T.plateau_kN(m));
            [yo, yn] = deal(zeros(size(F)));
            for i = 1:numel(F)
                mi = m & T.plateau_kN == F(i);
                switch q
                    case 1
                        yo(i) = max(T.main_comp_kN(mi & T.set == "offnominal"));
                        yn(i) = max(T.main_comp_kN(mi & T.set == "nominal"));
                    case 2
                        yo(i) = max(max(T.stroke_core_mm(mi & T.set == "offnominal"), T.stroke_pad_mm(mi & T.set == "offnominal"))) / 0.75;
                        yn(i) = max(max(T.stroke_core_mm(mi & T.set == "nominal"), T.stroke_pad_mm(mi & T.set == "nominal")));
                    case 3
                        yo(i) = max(T.acc_lp_g(mi & T.set == "offnominal"));
                        yn(i) = max(T.acc_lp_g(mi & T.set == "nominal"));
                end
            end
            h = plot(ax, F, yo, '-o', 'DisplayName', strrep(char(surf(s)), '_', ' '));
            plot(ax, F, yn, '--', 'Color', h.Color, 'HandleVisibility', 'off');
        end
        xlabel(ax, sprintf('%s core plateau force (kN)', place(p)));
        switch q
            case 1
                yline(ax, Pcr_main / opts.FS_buckle, 'r:', 'P_{cr}/2', 'HandleVisibility', 'off');
                ylabel(ax, 'peak main-strut compression (kN)');
            case 2, ylabel(ax, 'core length needed (mm); nominal crush dashed');
            case 3, ylabel(ax, 'peak CG accel, 100 Hz (g)');
        end
        legend(ax, 'Location', 'best');
    end
end
exportgraphics(f, fullfile(opts.OutDir, 'landing_sweep.png'), 'Resolution', 150);
close(f);
end

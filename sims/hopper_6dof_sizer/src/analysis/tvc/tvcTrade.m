function T = tvcTrade(opts)
%TVCTRADE  AG 9113 vs AG 9313 TVC trade and gimbal requirements.
%   T = tvcTrade()   runs both actuators through the sizer kinematics and
%   the PDR-slide-170 voltage-to-angle model, compares them against the
%   GNC gimbal requirements, and writes figures + tables to
%   outputs/figures/tvc and results/tvc_trade.mat.
%
%   If results/tvc_demand.mat exists (tvcDemandFlights) the 6DOF flight
%   demand is included.

arguments
    opts.OutDir (1,:) char = fullfile(hopperRootSafe(), 'outputs', 'figures', 'tvc')
    opts.DemandFile (1,:) char = fullfile(hopperRootSafe(), 'results', 'tvc_demand.mat')
end
if ~exist(opts.OutDir, 'dir'), mkdir(opts.OutDir); end
names = {'AG9113', 'AG9313'};
REQ = tvcRequirements();

%% ---- Arm-length study: rod end moved radially, tail stock moved so the
%      actuator sits at mid-stroke at 0 deg (L0 = datasheet M) ----
arms = (0.05:0.0025:0.16);
for a = 1:2
    for k = 1:numel(arms)
        z = solveRodEnd(names{a}, arms(k));
        P = tvcActuatorParams(names{a}, RodEndZ=z, FitLength=true);
        Aq(a, k) = quick(P, REQ); %#ok<AGROW>
    end
end

%% ---- Recommended arm: balance rate and force margin (both normalized) ----
for a = 1:2
    q = Aq(a, :);
    score = min([q.rate_load_dps] / REQ.rate_dps, [q.torque_margin] / REQ.force_margin);
    ok = score >= 1 & [q.res_deg] <= REQ.res_deg & [q.stroke_used] <= 0.9;
    score(~ok & any(ok)) = -inf;
    [~, k] = max(score);
    REC(a) = struct('arm_m', arms(k), 'feasible', ok(k), 'z', solveRodEnd(names{a}, arms(k))); %#ok<AGROW>
end

%% ---- Full evaluation: sizer layouts as drawn, and recommended layouts ----
cfg = {tvcActuatorParams('AG9113'), 'AG 9113, sizer layout'
       tvcActuatorParams('AG9313'), 'AG 9313, sizer layout'
       tvcActuatorParams('AG9113', RodEndZ=REC(1).z, FitLength=true), 'AG 9113, recommended'
       tvcActuatorParams('AG9313', RodEndZ=REC(2).z, FitLength=true), 'AG 9313, recommended'};
for c = 1:size(cfg, 1)
    E(c) = evaluate(cfg{c, 1}, REQ, cfg{c, 2}); %#ok<AGROW>
end

%% ---- Rotor-mass sensitivity (AG 9313, estimated m_r) ----
mrs = [50 100 200];
for k = 1:numel(mrs)
    P = tvcActuatorParams('AG9313');
    g = mrs(k) / P.m_r;  P.m_r = mrs(k);
    P.Kp = P.Kp*g; P.Ki = P.Ki*g; P.Kd = P.Kd*g;
    K = tvcKinematics(P, -8:0.5:8);
    s = stepMetrics(P, K, 0.5);  s8 = stepMetrics(P, K, 8);
    MR(k) = struct('m_r', mrs(k), 'rise05_ms', s.rise*1e3, 'os05', s.os, ...
                   'rise8_s', s8.rise, 'phase2Hz', sineResp(P, K, 0.3, 2).phase); %#ok<AGROW>
end

%% ---- Flight demand ----
DEM = [];
if exist(opts.DemandFile, 'file')
    L = load(opts.DemandFile);
    DEM = demandStats(L.D);
end

T = struct('REQ', REQ, 'E', E, 'arms', arms, 'Aq', Aq, 'MR', MR, 'DEM', DEM, 'REC', REC);
save(fullfile(hopperRootSafe(), 'results', 'tvc_trade.mat'), '-struct', 'T');
makeFigures(T, opts.OutDir);
printSummary(T);
end

%% =====================================================================
function E = evaluate(P, REQ, label)
E.label = label;
K  = tvcKinematics(P, -REQ.range_deg:0.1:REQ.range_deg);
Kf = tvcKinematics(P, -8:0.5:8);
E.name = P.name;  E.P = P;  E.TS_x_mm = P.TS_P(1)*1e3;  E.RE_z_mm = P.RE_P(3)*1e3;
E.L0_mm = P.L0*1e3;  E.M_mm = P.M*1e3;  E.fits = abs(P.L0 - P.M) < 0.005;
E.r_eff_mm = P.r_eff*1e3;
% stroke over the square +/- range envelope (both axes)
E.ext_mm = 1e3*max([K.dL_P(:); K.dL_Y(:)]);  E.ret_mm = -1e3*min([K.dL_P(:); K.dL_Y(:)]);
E.stroke_used = max(E.ext_mm, E.ret_mm) / (P.stroke*1e3);
% worst-case rod speed with both axes at the required rate (sizer method)
w = deg2rad(REQ.rate_dps);
vP = abs(K.dLdth_P)*w + abs(K.dLdps_P)*w;  vY = abs(K.dLdth_Y)*w + abs(K.dLdps_Y)*w;
E.v_req_mm_s = 1e3*max([vP(:); vY(:)]);
% load torque & force
E.tau_load = loadTorque(P, REQ);
r_min = min(abs(K.axis.r(abs(K.axis.delta) <= deg2rad(REQ.range_deg))));
E.F_load = E.tau_load / r_min;
E.torque_margin = P.Fmax / E.F_load;
E.tau_nom = loadTorque(P, REQ, P.e_nom);
E.torque_margin_nom = P.Fmax * r_min / E.tau_nom;
i_load = E.F_load / P.KF;
E.v_avail_mm_s = 1e3*(P.Vbus - P.R0*i_load) / P.Ke;
E.speed_margin = E.v_avail_mm_s / E.v_req_mm_s;
E.rate_load_dps = rad2deg((P.Vbus - P.R0*i_load)/P.Ke / P.r_eff);
E.rate_noload_dps = rad2deg(P.rate_noload);
E.res_deg = rad2deg(P.res_angle);
E.accel_cap_dps2 = rad2deg((P.KF*P.Imax - P.Fc - E.F_load) / (P.m_r + P.Jg/P.r_eff^2) / P.r_eff);
% dynamics
for A = [0.5 2 8]
    E.step.(sprintf('a%g', A*10)) = stepMetrics(P, Kf, A);
end
fr = [0.5 1 2 3 5 7 10];
for k = 1:numel(fr)
    s1 = sineResp(P, Kf, 0.3, fr(k));
    E.freq(k) = struct('f', fr(k), 'gain_dB', s1.gain_dB, 'phase', s1.phase); %#ok<AGROW>
end
s2 = sineResp(P, Kf, 0.3, 2);  E.phase2Hz = s2.phase;
% largest amplitude at the pitch crossover before rate saturation bites
Aset = [0.5 1 1.5 2 2.5 3 4];  fc = REQ.wc_rad_s/(2*pi);
ph0 = sineResp(P, Kf, 0.1, fc).phase;
E.A_lin_deg = NaN;
for A = Aset
    s = sineResp(P, Kf, A, fc);
    if s.phase < ph0 - 10 || s.gain_dB < -1, break; end
    E.A_lin_deg = A;
end
E.mass_pair_kg = 2*P.mass;
end

function q = quick(P, REQ)
K = tvcKinematics(P, -REQ.range_deg:0.5:REQ.range_deg);
q.r_eff = P.r_eff;
tau = loadTorque(P, REQ);
r_min = min(abs(K.axis.r(abs(K.axis.delta) <= deg2rad(REQ.range_deg))));
F = tau / r_min;  q.torque_margin = P.Fmax / F;
i = F / P.KF;
w = deg2rad(REQ.rate_dps);
vP = abs(K.dLdth_P)*w + abs(K.dLdps_P)*w;
q.rate_load_dps = REQ.rate_dps * ((P.Vbus - P.R0*i)/P.Ke) / max(vP(:));  % simultaneous-axes rate
q.stroke_used = 1e3*max(abs([K.dL_P(:); K.dL_Y(:)])) / (P.stroke*1e3);
q.res_deg = rad2deg(P.q_mech / P.r_eff);
end

function tau = loadTorque(P, REQ, e)
% gimbal hinge torque at the range edge, worst sign (stress offset default)
if nargin < 3, e = P.e_stress; end
d = deg2rad(REQ.range_deg);  a = deg2rad(REQ.accel_dps2);
tau = P.kh*d + P.tauf + P.thrust*e + P.Jg*a + P.me*REQ.axial_g*9.81*P.rcg*sin(d);
end

function m = stepMetrics(P, K, A)
S = tvcAxisSim(P, K, max(0.6, 0.15 + deg2rad(A)/P.rate_noload*1.6), @(t) deg2rad(A)*(t >= 0));
d = S.delta / deg2rad(A);
i10 = find(d >= 0.1, 1);  i90 = find(d >= 0.9, 1);
m.rise = S.t(i90) - S.t(i10);
m.os = 100*(max(d) - 1);
m.slew_dps = 0.8*A / m.rise;
m.peak_rate_dps = rad2deg(max(abs(S.rate)));
m.peak_accel_dps2 = rad2deg(max(abs(S.accel)));
m.peak_I = max(abs(S.i));
m.S = S;
end

function r = sineResp(P, K, A, f)
nc = max(3, ceil(f*1.0));  t_end = 0.2 + (nc + 2)/f;
S = tvcAxisSim(P, K, t_end, @(t) deg2rad(A)*sin(2*pi*f*t), decimate=20);
m = S.t > t_end - nc/f;
X = [sin(2*pi*f*S.t(m)), cos(2*pi*f*S.t(m))];
c = X \ S.delta(m);
r.gain_dB = 20*log10(hypot(c(1), c(2)) / deg2rad(A));
r.phase = rad2deg(atan2(c(2), c(1)));
end

function z = solveRodEnd(name, arm)
% rod-end radial offset that gives r_eff = arm with L0 = M
f = @(z) tvcActuatorParams(name, RodEndZ=z, FitLength=true).r_eff - arm;
z = fzero(f, [0.005 0.25]);
end

function D = demandStats(Dr)
D = struct('name', {}, 'max_deg', {}, 'max_rate_dps', {}, 'p99_rate_dps', {}, ...
           'max_accel_dps2', {}, 'p99_accel_dps2', {});
for k = 1:numel(Dr)
    if ~isempty(Dr(k).error) || isempty(Dr(k).t), continue; end
    t = Dr(k).t(:);  [t, iu] = unique(t);
    tt = (t(1):0.002:t(end))';
    ang = [interp1(t, Dr(k).pitch(iu), tt), interp1(t, Dr(k).yaw(iu), tt)];
    % the logged angle carries the 0.1 deg quantizer: smooth at ~10 Hz
    % before differentiating so steps of one count are not read as demand
    [b, a] = butter(2, 10/(500/2));
    ang_f = filtfilt(b, a, ang);
    rate = gradient(ang_f', 0.002)';  acc = gradient(rate', 0.002)';
    D(end+1) = struct('name', Dr(k).name, 'max_deg', rad2deg(max(abs(ang(:)))), ...
        'max_rate_dps', rad2deg(max(abs(rate(:)))), 'p99_rate_dps', rad2deg(prctile(abs(rate(:)), 99)), ...
        'max_accel_dps2', rad2deg(max(abs(acc(:)))), 'p99_accel_dps2', rad2deg(prctile(abs(acc(:)), 99))); %#ok<AGROW>
end
end

function makeFigures(T, outDir)
c = [0.13 0.27 0.53; 0.85 0.55 0.0; 0.13 0.27 0.53; 0.85 0.55 0.0];
ls = {'--', '--', '-', '-'};  nE = numel(T.E);  lab = {T.E.label};
% steps
f = figure('Visible', 'off', 'Position', [0 0 1200 420], 'Color', 'w', 'Theme', 'light');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
for A = [5 20 80]
    nexttile(tl); hold on; grid on
    for a = 1:nE
        s = T.E(a).step.(sprintf('a%g', A));
        plot(s.S.t*1e3, rad2deg(s.S.delta), ls{a}, 'Color', c(a,:), 'LineWidth', 1.5);
    end
    yline(A/10, ':'); xlabel('Time (ms)'); ylabel('Gimbal angle (deg)');
    title(sprintf('%.1f deg step', A/10)); legend(lab, 'Location', 'southeast');
end
exportgraphics(f, fullfile(outDir, 'tvc_steps.png'), 'Resolution', 200); close(f)
% frequency
f = figure('Visible', 'off', 'Position', [0 0 900 600], 'Color', 'w', 'Theme', 'light');
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
nexttile(tl); hold on; grid on; set(gca, 'XScale', 'log')
for a = 1:nE, semilogx([T.E(a).freq.f], [T.E(a).freq.gain_dB], ['o' ls{a}], 'Color', c(a,:), 'LineWidth', 1.5); end
ylabel('Gain (dB)'); title('0.3 deg sine response'); legend(lab)
nexttile(tl); hold on; grid on; set(gca, 'XScale', 'log')
for a = 1:nE, semilogx([T.E(a).freq.f], [T.E(a).freq.phase], ['o' ls{a}], 'Color', c(a,:), 'LineWidth', 1.5); end
xline(T.REQ.wc_rad_s/2/pi, ':', 'pitch crossover'); xlabel('Frequency (Hz)'); ylabel('Phase (deg)')
exportgraphics(f, fullfile(outDir, 'tvc_freq.png'), 'Resolution', 200); close(f)
% arm study
f = figure('Visible', 'off', 'Position', [0 0 1200 420], 'Color', 'w', 'Theme', 'light');
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
ax = T.arms*1e3;
nexttile(tl); hold on; grid on
for a = 1:2, plot(ax, [T.Aq(a,:).rate_load_dps], 'Color', c(a,:), 'LineWidth', 1.5); end
yline(T.REQ.rate_dps, 'k--', 'req'); xlabel('Moment arm r_{eff} (mm)'); ylabel('Rate at load, both axes (deg/s)')
legend({'AG 9113', 'AG 9313'}); for a = 1:2, xline(T.REC(a).arm_m*1e3, ':', 'Color', c(a,:), 'HandleVisibility', 'off'); end
nexttile(tl); hold on; grid on
for a = 1:2, plot(ax, [T.Aq(a,:).torque_margin], 'Color', c(a,:), 'LineWidth', 1.5); end
yline(T.REQ.force_margin, 'k--', 'req'); xlabel('Moment arm r_{eff} (mm)'); ylabel('Force margin (F_{rated}/F_{load})')
nexttile(tl); hold on; grid on
for a = 1:2, plot(ax, [T.Aq(a,:).res_deg], 'Color', c(a,:), 'LineWidth', 1.5); end
yline(T.REQ.res_deg, 'k--', 'req'); xlabel('Moment arm r_{eff} (mm)'); ylabel('Angle resolution (deg)')
exportgraphics(f, fullfile(outDir, 'tvc_arm_trade.png'), 'Resolution', 200); close(f)
end

function printSummary(T)
fprintf('\n==== TVC trade ====\n');
for a = 1:numel(T.E)
    E = T.E(a);
    fprintf('%s: RE z %.1f mm, TS x %.1f mm, L0 %.0f mm (mid %.0f, fits %d), r_eff %.1f mm, stroke used %.0f %% (+%.1f/-%.1f mm)\n', ...
        E.label, E.RE_z_mm, E.TS_x_mm, E.L0_mm, E.M_mm, E.fits, E.r_eff_mm, 100*E.stroke_used, E.ext_mm, E.ret_mm);
    fprintf('   torque stress %.1f / nom %.1f N m -> %.0f N, margin %.2fx (nom %.2fx); rod speed req %.1f vs %.1f mm/s (%.2fx)\n', ...
        E.tau_load, E.tau_nom, E.F_load, E.torque_margin, E.torque_margin_nom, E.v_req_mm_s, E.v_avail_mm_s, E.speed_margin);
    fprintf('   rate at load %.1f dps (no-load %.1f), accel cap %.0f dps2, res %.3f deg\n', ...
        E.rate_load_dps, E.rate_noload_dps, E.accel_cap_dps2, E.res_deg);
    fprintf('   0.5: rise %.0f ms OS %.1f%% pkacc %.0f | 2: %.0f ms | 8: %.0f ms slew %.1f pkrate %.1f pkacc %.0f | ph2Hz %.0f | A_lin %.1f\n', ...
        1e3*E.step.a5.rise, E.step.a5.os, E.step.a5.peak_accel_dps2, 1e3*E.step.a20.rise, 1e3*E.step.a80.rise, ...
        E.step.a80.slew_dps, E.step.a80.peak_rate_dps, E.step.a80.peak_accel_dps2, E.phase2Hz, E.A_lin_deg);
end
for a = 1:2
    fprintf('recommended arm %d: %.1f mm (feasible %d)\n', a, 1e3*T.REC(a).arm_m, T.REC(a).feasible);
end
disp(struct2table(T.MR));
if ~isempty(T.DEM), disp(struct2table(T.DEM)); end
end

function r = hopperRootSafe()
try
    r = hopperRoot();
catch
    r = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
end

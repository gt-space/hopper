function K = tvcKickTest(opts)
%TVCKICKTEST  Pitch-axis kick test: inject 2 Hz gimbal disturbance, check decay/sustain.
%
%   K = tvcKickTest()
%
%   Model: pitch double integrator + PD GNC + tvcEmUpdate (AG9113, 20.5 deg/s,
%   235 deg/s^2).  Same plant and GNC calibration as tvcDfMargins (PM_nominal=40 deg).
%   This is a pitch-axis ODE, not the full Simulink 6DOF — it is exactly the model
%   the DF analysis was computed for, so the kick results directly verify the DF
%   limit-cycle prediction (A_lc = 1.62 deg, |L_OL| = 0.995 at 2 Hz).
%
%   For each amplitude in A_kick_deg the sim runs:
%     [0, t_settle]           : closed-loop from rest (trim = theta=0)
%     [t_settle, t_kick_off]  : forced disturbance d(t) = A*sin(2pi*f*(t-t_settle))
%                               added to the GNC gimbal command
%     [t_kick_off, t_end]     : free response — watched for decay or sustain
%
%   Cases run with parfor. Each is a single nominal run, no MC.

arguments
    opts.A_kick_deg  (1,:) double = [1.0 1.5 2.0 2.5]   % deg
    opts.f_kick_Hz   (1,1) double = 2.0
    opts.N_kick      (1,1) double = 3        % forced cycles before release
    opts.t_settle    (1,1) double = 2.0      % s, let PD drive to trim
    opts.t_post      (1,1) double = 5.0      % s to observe after kick
    opts.PM_nominal  (1,1) double = 40       % deg, PD calibration
    opts.SaveFile    (1,:) char = fullfile(hopperRootSafe(),'results','tvc_kick_test.mat')
    opts.OutDir      (1,:) char = fullfile(hopperRootSafe(),'outputs','figures','tvc')
end
if ~exist(opts.OutDir,'dir'), mkdir(opts.OutDir); end
t0 = tic;

%% ---- Plant / GNC parameters ------------------------------------------
Iy = 41; Lca = 1.040 - 0.45; Th = 1116;
k_p = Th * Lca / Iy;        % 16.06 rad/s^2/rad
REQ = tvcRequirements();
wc = REQ.wc_rad_s; wn = REQ.wn_act; zeta = 0.7;

% PD gains (same calibration as tvcDfMargins)
PM0    = opts.PM_nominal;
Ta_wc  = wn^2 / ((1i*wc)^2 + 2*zeta*wn*(1i*wc) + wn^2);
ph_C   = PM0 - rad2deg(angle(Ta_wc));
Kd_Kp  = tand(ph_C) / wc;
C_norm = sqrt(Kd_Kp^2*wc^2 + 1);
Kp     = wc^2 / (k_p * C_norm * abs(Ta_wc));
Kd     = Kd_Kp * Kp;
fprintf('PD gains (PM_lin=%d deg): Kp=%.4f, Kd=%.4f\n', PM0, Kp, Kd);

%% ---- Actuator setup (AG9113, 20.5 deg/s, 235 deg/s^2) ----------------
P    = tvcActuatorParams('AG9113', Layout='recommended');
rlim = deg2rad(20.5);
alim = deg2rad(235);
T    = tvcEmParams(P, Offset=[P.e_nom; P.e_nom], Coulomb=true, Backlash=true, ...
         MaxAngle=deg2rad(20), AccelLimit=alim, RateLimit=rlim);
dt   = T.dt_c;   % 1 ms servo period

%% ---- Time grid -------------------------------------------------------
f          = opts.f_kick_Hz;
t_kick_on  = opts.t_settle;
t_kick_off = t_kick_on + opts.N_kick / f;
t_end      = t_kick_off + opts.t_post;
N          = round(t_end / dt);
t_vec      = (0:N)' * dt;

fprintf('[timing] setup: %.2f s\n', toc(t0));
fprintf('Running %d cases (parfor), each %.1f s @ %.0f Hz servo...\n', ...
    numel(opts.A_kick_deg), t_end, 1/dt);

%% ---- Parfor loop -----------------------------------------------------
n_cases     = numel(opts.A_kick_deg);
A_kick_arr  = opts.A_kick_deg;
results_arr = cell(1, n_cases);

t_sim = tic;
parfor ci = 1:n_cases
    A_kick = deg2rad(A_kick_arr(ci));

    St  = tvcEmInit();
    T_w = T;   % parfor slice

    theta_v    = zeros(N+1, 1);
    thetadot_v = zeros(N+1, 1);
    delta_v    = zeros(N+1, 1);
    cmd_v      = zeros(N+1, 1);

    for k = 1:N
        tk  = t_vec(k);
        th  = theta_v(k);
        thd = thetadot_v(k);

        % Stabilising PD (drives theta → 0)
        gnc = -(Kp * th + Kd * thd);

        % Additive kick disturbance (sine burst at 2 Hz)
        if tk >= t_kick_on && tk < t_kick_off
            gnc = gnc + A_kick * sin(2*pi*f * (tk - t_kick_on));
        end
        cmd_v(k) = gnc;

        % Actuator step (pitch axis; yaw held at 0)
        [St, dg, ~, ~] = tvcEmUpdate(St, [gnc; 0], Th, T_w);
        delta_v(k+1) = dg(1);

        % Rigid body (semi-implicit Euler)
        alpha          = k_p * dg(1);
        thetadot_v(k+1) = thd + dt * alpha;
        theta_v(k+1)    = th  + dt * thetadot_v(k+1);
    end
    cmd_v(N+1) = -(Kp*theta_v(N+1) + Kd*thetadot_v(N+1));

    % ---- Assess decay / sustain in post-kick window -------------------
    %   Compare RMS in first 1 s and last 1 s of post-kick window
    %   (skip a 0.5 s transient after kick ends)
    t_obs_start = t_kick_off + 0.5;
    m_early = t_vec >= t_obs_start         & t_vec < t_obs_start + 1.0;
    m_late  = t_vec >= t_end - 1.0;

    rms_early_deg = rad2deg(rms(theta_v(m_early)));
    rms_late_deg  = rad2deg(rms(theta_v(m_late)));

    results_arr{ci} = struct( ...
        'A_kick_deg',    A_kick_arr(ci), ...
        't',             t_vec, ...
        'theta_deg',     rad2deg(theta_v), ...
        'thetadot_dps',  rad2deg(thetadot_v), ...
        'delta_deg',     rad2deg(delta_v), ...
        'cmd_deg',       rad2deg(cmd_v), ...
        'rms_early_deg', rms_early_deg, ...
        'rms_late_deg',  rms_late_deg);
end
fprintf('[timing] parfor sim: %.2f s\n', toc(t_sim));

%% ---- Report ----------------------------------------------------------
fprintf('\n=== Kick test: AG9113 20.5°/s 235°/s², f=%.1f Hz, PM_lin=%d° ===\n', f, PM0);
fprintf('  DF prediction: A_lc=1.62°, |L_OL|=0.995 (marginal LC at 2 Hz)\n\n');
fprintf('%-11s  %10s  %9s  %s\n', 'A_kick(°)', 'RMS_early°', 'RMS_late°', 'Verdict');
for ci = 1:n_cases
    r = results_arr{ci};
    ratio = r.rms_late_deg / max(r.rms_early_deg, 1e-6);
    if ratio < 0.5
        verdict = 'DECAYS';
    elseif ratio < 1.5
        verdict = 'SUSTAINS (marginal)';
    else
        verdict = 'GROWS';
    end
    fprintf('%-11.2f  %10.3f  %9.3f  %s (ratio=%.2f)\n', ...
        r.A_kick_deg, r.rms_early_deg, r.rms_late_deg, verdict, ratio);
end

%% ---- Figure: 2×2 pitch time history ----------------------------------
t_fig = tic;
clrs = [0 0.45 0.74; 0.85 0.33 0.1; 0.47 0.67 0.19; 0.49 0.18 0.56];

fig = figure('Visible','off','Position',[50 50 1100 720],'Color','w','Theme','light');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('Kick test: AG9113 20.5°/s 235°/s²,  f=%.0f Hz,  PM_{lin}=%d°', f, PM0));

for ci = 1:n_cases
    r   = results_arr{ci};
    ax  = nexttile(tl);
    hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    xline(ax, t_kick_on,  'Color',[0.5 0.5 0.5], 'LineStyle',':', 'Label','kick on',  'LabelVerticalAlignment','bottom');
    xline(ax, t_kick_off, 'Color',[0.5 0.5 0.5], 'LineStyle','--','Label','kick off', 'LabelVerticalAlignment','bottom');
    plot(ax, r.t, r.theta_deg, 'Color', clrs(ci,:), 'LineWidth', 1.5);
    xlabel(ax,'Time (s)');  ylabel(ax,'Pitch angle (deg)');
    title(ax, sprintf('A_{kick}=%.1f°  (early %.3f° / late %.3f°)', ...
        r.A_kick_deg, r.rms_early_deg, r.rms_late_deg));
end

exportgraphics(fig, fullfile(opts.OutDir,'tvc_kick_test.png'), 'Resolution',200);
fprintf('[timing] figure exported: %.2f s\n', toc(t_fig));
K.fig = fig;

%% ---- Save (data only, no figure handles) -----------------------------
t_sv = tic;
K.results      = results_arr;
K.Kp           = Kp;  K.Kd = Kd;  K.PM_nominal = PM0;
K.A_kick_deg   = opts.A_kick_deg;  K.f_kick_Hz = f;
K.t_kick_on    = t_kick_on;  K.t_kick_off = t_kick_off;
K.t_end        = t_end;
K_save         = rmfield(K, 'fig');
save(opts.SaveFile, '-struct', 'K_save', '-v7.3');
fprintf('[timing] save: %.2f s\n', toc(t_sv));
fprintf('[timing] TOTAL: %.2f s\n', toc(t0));
end

function r = hopperRootSafe()
try, r = hopperRoot();
catch, r = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))); end
end

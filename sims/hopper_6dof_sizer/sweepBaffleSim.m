function BS = sweepBaffleSim(overrides)
% Closed-loop baffle sizing: sweep baffle count / width per tank in
% hopper_6dof_NED_v2 and pick the fewest baffles that keep altitude and
% pitch within limits, with margin on the damping model.
%
% Run sim_setup first. This is a function so nothing it creates can
% overwrite the model's base workspace variables (e.g. g); it only reads
% from base and runs generateBaffleLUT there.
%
% Margin: every design must pass with all damping tables scaled by each
% of BS.zeta_scales. Runs are not bit-repeatable and a sparse design can
% diverge on only some runs, so each case is run BS.n_repeat times and
% scored on its worst run. 2 repeats proved too few: designs that passed
% 6/6 runs diverged in a 5-repeat check at 0.5x damping.
% One tank is swept at a time with the other held at its analytic
% (generateBaffleLUT) design; the chosen pair is then re-checked together.
%
% Raw metrics are saved to baffle_sim_sweep.mat as soon as the sweep
% finishes. To re-score with different limits without re-simulating:
%   load baffle_sim_sweep.mat; BS.pitch_lim = 6; BS = scoreBaffleSweep(BS);

IN           = evalin('base', 'IN');
t            = evalin('base', 't');
z_trajectory = evalin('base', 'z_trajectory');

BS.mdl           = 'hopper_6dof_NED_v2';
BS.pitch_lim     = 5;      % max |pitch - 90| above BS.alt_split, deg
BS.pitch_lim_low = 10;     % same, final descent below BS.alt_split (Inf to ignore)
BS.alt_split     = 1;      % m
BS.alt_lim       = 1;      % max |altitude - reference|, m
BS.max_alt_min   = IN.mission.target_altitude;
BS.zeta_scales   = [1 0.5 2];
BS.n_repeat      = 5;
BS.Nb_vec        = 2:2:36;
BS.w_frac_vec    = [0.3 0.4 0.5];
BS.t             = 0.5e-3; % thickness only reduces damping; keep minimum
BS.n_workers     = 8;

% optional overrides, e.g. sweepBaffleSim(struct('Nb_vec', 10:2:20))
if nargin > 0
    for f = fieldnames(overrides)'
        BS.(f{1}) = overrides.(f{1});
    end
end

% analytic design (also sets tank lengths, densities, viscosities, LUTs)
evalin('base', ['baffle_mode = ''sweep''; baffle_export = false; generateBaffleLUT; ' ...
                'close all; clear baffle_mode baffle_export']);
baffle_opts     = evalin('base', 'baffle_opts');
tank_r          = evalin('base', 'tank_r');
ox_mass_profile = evalin('base', 'ox_mass_profile');
fu_mass_profile = evalin('base', 'fu_mass_profile');
ox_baffle       = evalin('base', 'ox_baffle');
fu_baffle       = evalin('base', 'fu_baffle');
BS.slosh_amp = baffle_opts.slosh_amp;
BS.tank_r    = tank_r;
BS.tank(1) = struct('name', 'ox', 'rho', evalin('base', 'ox_density'), 'nu', evalin('base', 'ox_nu'), ...
    'h', evalin('base', 'ox_tank_h'), 'm', ox_mass_profile, 'analytic', ox_baffle);
BS.tank(2) = struct('name', 'fu', 'rho', evalin('base', 'fu_density'), 'nu', evalin('base', 'fu_nu'), ...
    'h', evalin('base', 'fu_tank_h'), 'm', fu_mass_profile, 'analytic', fu_baffle);

BS.t_ref   = t(:);
BS.alt_ref = -z_trajectory(:);   % z_trajectory is NED

% ---- candidate list: stage 1 sweeps ox, stage 2 sweeps fu -------------
BS.cand = struct('stage', {}, 'design', {});
for bs_s = 1:2
    for bs_n = BS.Nb_vec
        for bs_w = BS.w_frac_vec
            d = [BS.tank.analytic];
            d = struct('Nb', {d.Nb}, 'w', {d.w}, 't', {d.t});
            d(bs_s) = struct('Nb', bs_n, 'w', bs_w * tank_r, 't', BS.t);
            BS.cand(end+1) = struct('stage', bs_s, 'design', d);
        end
    end
end

BS.metrics = runBaffleCases(BS, BS.cand);
save('baffle_sim_sweep.mat', 'BS');

% ---- pick fewest passing baffles per tank, then confirm together -------
BS = scoreBaffleSweep(BS);

BS.final = runBaffleCases(BS, struct('stage', 3, 'design', BS.final_design));
BS = scoreBaffleSweep(BS);
save('baffle_sim_sweep.mat', 'BS');

fprintf('\n=== Closed-loop baffle sizing (pitch <= %.1f deg above %.1f m, <= %.1f deg below, alt err <= %.2f m, damping x%s, %d repeats) ===\n', ...
    BS.pitch_lim, BS.alt_split, BS.pitch_lim_low, BS.alt_lim, mat2str(BS.zeta_scales), BS.n_repeat);
for bs_s = 1:2
    d = BS.final_design(bs_s); a = BS.tank(bs_s).analytic;
    fprintf('%s: Nb = %d, w = %.1f mm, t = %.1f mm   (analytic design: Nb = %d, w = %.1f mm)  stage pass: %d\n', ...
        BS.tank(bs_s).name, d.Nb, 1e3 * d.w, 1e3 * d.t, a.Nb, 1e3 * a.w, BS.stage_pass(bs_s));
end
fprintf('combined: worst pitch err = %.3f deg (high) / %.3f deg (low), worst alt err = %.3f m, min peak alt = %.2f m -> %s\n', ...
    max(BS.final.pitch_err_high(:)), max(BS.final.pitch_err_low(:)), max(BS.final.alt_err(:)), ...
    min(BS.final.max_alt(:)), string(BS.final.pass));

% ---- export LUT for the chosen design (only if it passed together) ------
if ~BS.final.pass
    warning('sweepBaffleSim:notExported', ...
        'Combined design failed its check; baffle_lut.csv left unchanged (analytic design).');
    return
end
[ox_lateral_damping_ratios, ox_axial_damping_ratios] = baffleLUT(BS.tank(1), BS.final_design(1), BS.tank_r, BS.slosh_amp);
[fu_lateral_damping_ratios, fu_axial_damping_ratios] = baffleLUT(BS.tank(2), BS.final_design(2), BS.tank_r, BS.slosh_amp);

baffle_lut = table(ox_mass_profile(:), ox_lateral_damping_ratios(:), ox_axial_damping_ratios(:), ...
    fu_mass_profile(:), fu_lateral_damping_ratios(:), fu_axial_damping_ratios(:), ...
    'VariableNames', {'ox_mass_profile', 'ox_lateral_damping_ratios', 'ox_axial_damping_ratios', ...
                      'fu_mass_profile', 'fu_lateral_damping_ratios', 'fu_axial_damping_ratios'});
writetable(baffle_lut, 'baffle_lut.csv');
baffle_sim_sweep = BS;
ox_tank_h = BS.tank(1).h; fu_tank_h = BS.tank(2).h;
save('baffle_lut.mat', 'baffle_lut', 'baffle_sim_sweep', 'baffle_opts', 'tank_r', 'ox_tank_h', 'fu_tank_h');
fprintf('Baffle LUT for the closed-loop design written to baffle_lut.csv / baffle_lut.mat\n');
end


function M = runBaffleCases(BS, cand)
% Run every candidate x damping scale x repeat in parallel.
nS = numel(BS.zeta_scales); nR = BS.n_repeat; nC = numel(cand);
in = repmat(Simulink.SimulationInput(BS.mdl), nC * nS * nR, 1);
t_ref = BS.t_ref; alt_ref = BS.alt_ref; alt_split = BS.alt_split;
k = 0;
for c = 1:nC
    [ol, oa] = baffleLUT(BS.tank(1), cand(c).design(1), BS.tank_r, BS.slosh_amp);
    [fl, fa] = baffleLUT(BS.tank(2), cand(c).design(2), BS.tank_r, BS.slosh_amp);
    for s = 1:nS
        g = BS.zeta_scales(s);
        for r = 1:nR
            k = k + 1;
            in(k) = in(k).setVariable('ox_lateral_damping_ratios', g * ol) ...
                         .setVariable('ox_axial_damping_ratios',   g * oa) ...
                         .setVariable('fu_lateral_damping_ratios', g * fl) ...
                         .setVariable('fu_axial_damping_ratios',   g * fa) ...
                         .setPostSimFcn(@(so) struct('m', baffleSimMetrics(so, t_ref, alt_ref, alt_split)));
        end
    end
end

if isempty(gcp('nocreate')), parpool('Processes', BS.n_workers); end
out = parsim(in, 'TransferBaseWorkspaceVariables', 'on', 'UseFastRestart', 'on', ...
    'ShowProgress', 'on');

fields = {'pitch_err', 'pitch_err_high', 'pitch_err_low', 'alt_err', 'max_alt'};
k = 0;
for c = 1:nC
    for f = fields, M(c).(f{1}) = NaN(nS, nR); end %#ok<AGROW>
    M(c).ok = false(nS, nR);
    for s = 1:nS
        for r = 1:nR
            k = k + 1;
            if isempty(out(k).ErrorMessage) && out(k).m.ok
                for f = fields, M(c).(f{1})(s, r) = out(k).m.(f{1}); end
                M(c).ok(s, r) = true;
            end
        end
    end
end
end

function [lat, ax] = baffleLUT(tank, d, tank_r, slosh_amp)
lat = lateralDampingCalc(tank.m, tank.rho, tank.nu, tank_r, d.w, d.Nb, tank.h, d.t, slosh_amp);
ax  = axialDampingCalc(tank.m, tank.rho, tank.nu, tank_r, d.w, d.Nb, tank.h, d.t, slosh_amp);
end

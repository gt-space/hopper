function C = tvcClampCampaign(opts)
%TVCCLAMPCAMPAIGN  6DOF TVC trade with the drive's rate/accel clamps applied.
%   C = tvcClampCampaign() flies models/hopper_6dof.slx (TVC_MODEL variant,
%   electromechanical actuator on the recommended layout) for the nominal
%   scenario plus the three stress cases of tvcLoopCampaign (12 m/s wind,
%   2.5 cm thrust-line offset, both), for each clamp config below, in
%   parallel. No dispersed seeds. Saves results/tvc_clamp_mc.mat.
%
%     AG 9113, no clamps                 baseline: motor/load limits only
%     AG 9113, 20.5 deg/s, 235 deg/s^2   TVC team's numbers at +/-8 deg
%     AG 9113, 20.5 deg/s, 450 deg/s^2   GNC desired upper end
%     AG 9313, 13.5 deg/s, 640 deg/s^2   TVC team's numbers at +/-8 deg
%
%   C = tvcClampCampaign(Configs=cfg) runs your own; each config is a
%   struct with name, model ("ag9113"/"ag9313"), layout, offset ([] =
%   follow the scenario, 1 cm floor), rate_dps, accel_dps2.
%   C = tvcClampCampaign(Seeds=50) adds 50 dispersed seeds.
%
%   Per run (C.runs): landing offset/requirements, peak gimbal rate and
%   accel vs clamp (rate_over_lim, accel_over_lim), % flight time at each
%   clamp, tracking error, tilt (max_tilt, rms_tilt_deg), and the sustained
%   oscillation check (osc_sustained, osc_amp_deg, osc_freq_hz, osc_dur_s).

arguments
    opts.Configs struct = clampConfigs()
    opts.Seeds (1,1) double = 0
    opts.OutputFile (1,:) char = hopperRoot('results', 'tvc_clamp_mc.mat')
    opts.Workers (1,1) double = NaN
end
C = tvcLoopCampaign(opts.Seeds, Configs=opts.Configs, OutputFile=opts.OutputFile, ...
                    Workers=opts.Workers, Stress=true);
cols = {'case', 'config', 'success', 'meets_requirements', 'landing_offset', 'max_tilt', 'rms_tilt_deg', ...
    'max_rate_dps', 'rate_lim_dps', 'rate_over_lim', 'at_rate_lim_pct', ...
    'max_accel_dps2', 'accel_lim_dps2', 'accel_over_lim', 'at_accel_lim_pct', ...
    'rms_err_deg', 'max_err_deg', 'osc_sustained', 'osc_amp_deg', 'osc_freq_hz', 'osc_dur_s', 'error'};
disp(sortrows(C.runs(:, cols), {'case', 'config'}));
end

function cfg = clampConfigs()
cfg = struct('name', {}, 'model', {}, 'layout', {}, 'offset', {}, 'rate_dps', {}, 'accel_dps2', {});
cfg(end+1) = struct('name', "AG 9113, no clamps (baseline)", 'model', "ag9113", 'layout', "recommended", 'offset', [], 'rate_dps', Inf, 'accel_dps2', Inf);
cfg(end+1) = struct('name', "AG 9113, 20.5 dps, 235 dps2", 'model', "ag9113", 'layout', "recommended", 'offset', [], 'rate_dps', 20.5, 'accel_dps2', 235);
cfg(end+1) = struct('name', "AG 9113, 20.5 dps, 450 dps2", 'model', "ag9113", 'layout', "recommended", 'offset', [], 'rate_dps', 20.5, 'accel_dps2', 450);
cfg(end+1) = struct('name', "AG 9313, 13.5 dps, 640 dps2", 'model', "ag9313", 'layout', "recommended", 'offset', [], 'rate_dps', 13.5, 'accel_dps2', 640);
end

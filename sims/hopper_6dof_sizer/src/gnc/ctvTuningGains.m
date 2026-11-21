function G = ctvTuningGains(knobs)
%CTVTUNINGGAINS  Gain tables for the CONTROLLER=lqr_ctv_tuning variant.
%   G = ctvTuningGains()          the tuned tables (the defaults below)
%   G = ctvTuningGains(knobs)     other weights, a struct with any of
%
%   The defaults hold the CTV's 0.5 m landing offset on the 5 seeded Monte
%   Carlo scenarios (inputs/mc_params_ctv.json, seed 42) flying on the GMEKF
%   estimate (NAVIGATION=gmekf): 0.09-0.34 m, against 0.3-10 m with the
%   nominal weights, of which two runs crashed. Measured 2026-11-20.
%       scale        pitch/yaw TVC effort weight multiplier   (default 30)
%                    raising it lowers the attitude-loop gain, which stops
%                    the command fighting TVC lag, quantization and rate
%       rateScale    body-rate weight multiplier              (default 10)
%                    more weight on the measured rates is more damping
%       lateralScale x, y and vx, vy weight multiplier        (default 1000)
%                    more weight is a stiffer station-keeping loop, which
%                    holds the landing offset down in wind
%       thrustScale  thrust effort weight multiplier          (default 1)
%       vehicle      design for this vehicle                  (default the
%                    current hopperVehicle)
%
%   Designs the time-scheduled LQR gains along the reference as sim_setup
%   does, with only those weight changes; the CTV's roll sizing, roll angle
%   weight and thrust scaling stay as inputs/ctv_inputs.m sets them.
%
%   Returns K1flat_ctv and K2grid_ctv, which subsystems/gnc/LQR_CTV_Tuning.slx
%   reads. Cached per vehicle and weights for the session, since every
%   scenario of a Monte Carlo shares the nominal vehicle's controller.

D = struct('scale', 30, 'rateScale', 10, 'lateralScale', 1000, 'thrustScale', 1, 'vehicle', "");
if nargin < 1 || isempty(knobs), knobs = struct(); end
for f = fieldnames(D)'
    if ~isfield(knobs, f{1}) || isempty(knobs.(f{1})), knobs.(f{1}) = D.(f{1}); end
end
if strlength(string(knobs.vehicle)) == 0, knobs.vehicle = hopperVehicle(); end

persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%s|%.6g|%.6g|%.6g|%.6g', knobs.vehicle, knobs.scale, knobs.rateScale, ...
    knobs.lateralScale, knobs.thrustScale);
if isKey(cache, key)
    G = cache(key);
    return
end

saved = evalin('base', 'who');
savedVals = cellfun(@(n) evalin('base', n), saved, 'uni', 0);
restore = onCleanup(@() restoreBase(saved, savedVals));
oldVehicle = hopperVehicle(knobs.vehicle);
restoreVehicle = onCleanup(@() hopperVehicle(oldVehicle));
evalin('base', 'clearvars');
evalin('base', 'sim_setup');

IN = evalin('base', 'IN');
if isfield(IN, 'gnc') && isfield(IN.gnc, 'lqr')
    L = IN.gnc.lqr;
    L.R_diag(1)   = L.R_diag(1) * knobs.thrustScale;
    L.R_diag(2:3) = L.R_diag(2:3) * knobs.scale;
    L.rate_weight_scale = knobs.rateScale;
    L.lateral_weight_scale = knobs.lateralScale;
    IN.gnc.lqr = L;
    assignin('base', 'IN', IN);
    % LinerizationMaster is a script and grows Fref, so a second run in the
    % same workspace appends to the last one's
    evalin('base', 'clear Fref');
    evalin('base', 'LinerizationMaster');
    G.K1flat_ctv = evalin('base', 'K1flat');
    G.K2grid_ctv = evalin('base', 'K2grid');
else
    % a vehicle with no CTV tuning (the hopper): the nominal tables, so the
    % choice flies the same controller instead of failing
    G.K1flat_ctv = evalin('base', 'K1flat');
    G.K2grid_ctv = evalin('base', 'K2grid');
end
cache(key) = G;
end

function restoreBase(names, vals)
evalin('base', 'clearvars');
for k = 1:numel(names)
    assignin('base', names{k}, vals{k});
end
end

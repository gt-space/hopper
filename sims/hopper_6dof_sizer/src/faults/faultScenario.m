function faults = faultScenario(varargin)
%FAULTSCENARIO  Fault table for the FAULT_INJECTION variant.
%   faults = faultScenario()                     no faults
%   faults = faultScenario(f1, f2, ...)          each f a cell:
%       {target, mode, t_on}
%       {target, mode, t_on, t_off}
%       {target, mode, t_on, t_off, magnitude}
%       {target, mode, t_on, t_off, magnitude, param}
%
%   target  'accel' 'gyro' 'baro' 'lidar' 'gps' 'mag'
%           'gimbal_pitch' 'gimbal_yaw' 'throttle' 'rcs'
%   mode    'loss' 'bias' 'drift' 'saturation' 'stuck' 'noise' 'stale'
%           'poor_fix' 'slow' 'hardover' 'dead' 'reduced'
%   t_off   default inf (permanent)
%   param   axis (accel/gyro), beam (lidar), satellites (poor_fix)
%
%   faults = faultScenario({'lidar', 'loss', 20, inf, 0, 2}, {'gimbal_pitch', 'stuck', 8})
%   out = sim(buildSimInput('hopper_6dof', [], W, variantConfig('FAULT_INJECTION', 'on', 'FDIR', 'on')))
%   with W.faults = faults (see docs/FDIR_ARCHITECTURE.md).
%
% The table is plain numbers ([target mode t_on t_off magnitude param]
% per row, fault_codes), so a Monte Carlo can disperse onsets and sizes.
% An empty scenario still has one inert row: Simulink needs a fixed size.

F = fault_codes();
targets = struct('accel', F.T_ACCEL, 'gyro', F.T_GYRO, 'baro', F.T_BARO, 'lidar', F.T_LIDAR, ...
    'gps', F.T_GPS, 'mag', F.T_MAG, 'gimbal_pitch', F.T_GIMBAL_P, 'gimbal_yaw', F.T_GIMBAL_Y, ...
    'throttle', F.T_THROTTLE, 'rcs', F.T_RCS);
modes = struct('loss', F.M_LOSS, 'bias', F.M_BIAS, 'drift', F.M_DRIFT, 'saturation', F.M_SATURATION, ...
    'stuck', F.M_STUCK, 'noise', F.M_NOISE, 'stale', F.M_STALE, 'poor_fix', F.M_POOR_FIX, ...
    'slow', F.M_SLOW, 'hardover', F.M_HARDOVER, 'dead', F.M_DEAD, 'reduced', F.M_REDUCED);

faults = zeros(max(1, nargin), F.NCOL);
faults(:, 3) = inf;   % the inert row never turns on
faults(:, 4) = inf;
for k = 1:nargin
    f = varargin{k};
    if ~iscell(f) || numel(f) < 3
        error('faultScenario:row', 'Fault %d: give {target, mode, t_on, ...}.', k);
    end
    tg = lower(char(f{1})); md = lower(char(f{2}));
    if ~isfield(targets, tg)
        error('faultScenario:target', 'Unknown target ''%s''. Targets: %s.', tg, strjoin(fieldnames(targets)', ', '));
    end
    if ~isfield(modes, md)
        error('faultScenario:mode', 'Unknown mode ''%s''. Modes: %s.', md, strjoin(fieldnames(modes)', ', '));
    end
    row = [targets.(tg), modes.(md), f{3}, inf, 0, 0];
    for j = 4:min(6, numel(f))
        row(j) = f{j};
    end
    faults(k, :) = row;
end
end

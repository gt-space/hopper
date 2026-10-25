function [raw_fail, monitor_valid, meas_bus, mon_state] = fdir_sensor_monitor( ...
        sensor_data, ekf_ref, source_mask, ground_reset, scfg, mon_state)
% FDIR_SENSOR_MONITOR  One sensor, one frame, four checks.
%
% Shared by all four sensor monitors; the per-sensor thresholds come in
% through scfg (see fdir_config).
%
%   raw_fail      4x1 logical, in order [loss; drift; saturation; stuck]
%   monitor_valid false when the source is masked out
%   meas_bus      the measurement, or NaN when the monitor is not valid
%   mon_state     stuck-check state, feed back in next frame
%                 (see fdir_monitor_state)
%
%#codegen

n = numel(sensor_data);
x = sensor_data(:);
ref = ekf_ref(:);

raw_fail = false(4, 1);

% A masked source produces no measurement and no fault reports. Every output
% is assigned on this path as well as the normal one.
if source_mask && ~ground_reset
    monitor_valid = false;
    meas_bus = NaN(n, 1);
    mon_state.initialized = false;
    mon_state.repeat_count = 0;
    return;
end

monitor_valid = true;
meas_bus = x;

% Check 1: loss of data. NaN and Inf only - a real zero is a valid reading
% (z = 0 sitting on the pad), so it must not count as loss.
raw_fail(1) = any(isnan(x) | isinf(x));

% Check 2: drift or bias against the EKF reference, per channel.
% Short-circuit when the data is already junk, so one NaN does not also
% raise drift and saturation.
if ~raw_fail(1)
    raw_fail(2) = any(abs(x - ref) > scfg.drift(:));

    % Check 3: saturation. Out of range is low OR high on any channel.
    raw_fail(3) = any(x < scfg.range_lo(:) | x > scfg.range_hi(:));
end

% Check 4: stuck sensor - bit-identical reading for stuck_frames frames in a
% row. The run length counts the current sample, so a fresh monitor starts
% at 1 and a sensor that never updates trips on frame stuck_frames.
if ~mon_state.initialized
    mon_state.repeat_count = 1;
    mon_state.initialized = true;
elseif isequal(x, mon_state.last_value(:))
    mon_state.repeat_count = mon_state.repeat_count + 1;
else
    mon_state.repeat_count = 1;
end
mon_state.last_value = x;

raw_fail(4) = mon_state.repeat_count >= scfg.stuck_frames;
end

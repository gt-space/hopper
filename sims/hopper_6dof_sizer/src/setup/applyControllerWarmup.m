function W = applyControllerWarmup(W, warmup_time)
%APPLYCONTROLLERWARMUP Delay the reference schedule and zero actuator command.
if ~isnumeric(warmup_time) || ~isscalar(warmup_time) || ~isfinite(warmup_time) || warmup_time < 0
    error('applyControllerWarmup:invalidTime', 'Warm-up time must be a finite, nonnegative scalar.');
end
if ~isfield(W, 'tgrid') || ~isfield(W, 'unom') || isempty(W.tgrid) || isempty(W.unom)
    error('applyControllerWarmup:missingSchedule', 'Workspace must contain a nonempty tgrid and unom schedule.');
end

if isfield(W, 'controller_start_time')
    previous_time = W.controller_start_time;
else
    previous_time = 0;
end
time_shift = warmup_time - previous_time;
if time_shift ~= 0
    W.tgrid = W.tgrid + time_shift;
    series_names = {'T_ts', 'T_ref_ts', 'z_ts', 'z_ref_ts', 'vz_ts', 'vz_ref_ts'};
    for k = 1:numel(series_names)
        name = series_names{k};
        if isfield(W, name) && isa(W.(name), 'timeseries')
            ts = W.(name);
            ts.Time = ts.Time + time_shift;
            W.(name) = ts;
        end
    end
    if isfield(W, 'scenario') && isa(W.scenario, 'Simulink.SimulationData.Dataset')
        for k = 1:W.scenario.numElements
            ts = W.scenario.getElement(k);
            if isa(ts, 'timeseries')
                ts.Time = ts.Time + time_shift;
                W.scenario = W.scenario.setElement(k, ts);
            end
        end
    end
end

W.controller_start_time = warmup_time;
W.u_hold_command = zeros(size(W.unom, 2), 1);
end

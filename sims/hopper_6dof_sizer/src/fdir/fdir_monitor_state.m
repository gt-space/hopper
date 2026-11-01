function mon_state = fdir_monitor_state(n)
% FDIR_MONITOR_STATE  Fresh stuck-check state for an n-channel sensor.
%
%#codegen

mon_state = struct( ...
    'last_value',   zeros(n, 1), ...
    'repeat_count', 0, ...
    'initialized',  false);
end

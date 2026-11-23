function W = ctrlRateWorkspace(W, period, latency)
%CTRLRATEWORKSPACE  Controller timing for CONTROLLER=lqr_ctv_tuning.
%   W = ctrlRateWorkspace(W)                  every solver step (the default;
%                                             keeps values already in W)
%   W = ctrlRateWorkspace(W, period, latency) the flight computer's timing:
%                                             the scheduler reads the state
%                                             and clock at each period's
%                                             start and its command lands
%                                             latency later, held for the
%                                             rest of the period
%
%   Sets the sample times of the holds addCtvTuningController puts around
%   the scheduler in subsystems/gnc/LQR_CTV_Tuning.slx:
%     ctrl_ts_in   state and clock hold: -1 (inherited) or [period, 0]
%     ctrl_ts_out  command hold: -1 or [period, latency]
%   latency must be at least one solver step (1 ms) and under the period.
%
%   Example, luna flight2's 100 Hz loop with a 1 ms command latency:
%     V = variantConfig('VEHICLE', 'ctv');
%     W = ctrlRateWorkspace(hopperWorkspace([], V), 0.01, 1e-3);
%     out = sim(buildSimInput('hopper_6dof', [], W, V));

if nargin < 2 || isempty(period)
    if ~isfield(W, 'ctrl_ts_in'), W.ctrl_ts_in = -1; end
    if ~isfield(W, 'ctrl_ts_out'), W.ctrl_ts_out = -1; end
    return
end
assert(period > 0 && latency > 0 && latency < period, 'ctrlRateWorkspace:timing', ...
    'need 0 < latency < period (got period %g, latency %g)', period, latency);
W.ctrl_ts_in = [period, 0];
W.ctrl_ts_out = [period, latency];
end

function W = fswWorkspace(W)
%FSWWORKSPACE  Variables for CONTROLLER=fsw (flight software flies).
%   W = fswWorkspace(W) adds what the fsw choice reads, unless W already
%   has it:
%     fsw_period  s, the flight software's control period: the root inport
%                 u_fsw holds each command for one period (100 Hz,
%                 luna flight2)
%     fsw_latency s, when in each period the command takes effect: the
%                 flight software reads the state at the period's start
%                 and its command lands one plant step (1 ms) later
%     fsw_u       the command replayed on u_fsw in a simulation, rows of
%                 [t, thrust N, TVC pitch rad, TVC yaw rad, RCS torque N m].
%                 Held between rows, not interpolated. Default: the nominal
%                 feedforward unom along tgrid, open loop. Set it to a
%                 flight-software log to replay that run through Simulink
%                 (fsw_bridge/verify_plant_c).
%
%   Generated code (fsw_bridge/build_plant_c) ignores fsw_u: the flight
%   software writes u_fsw every period.

if ~isfield(W, 'fsw_period')
    W.fsw_period = 0.01;
end
if ~isfield(W, 'fsw_latency')
    W.fsw_latency = 1e-3;
end
if ~isfield(W, 'fsw_u')
    W.fsw_u = [W.tgrid(:), W.unom];
end
end

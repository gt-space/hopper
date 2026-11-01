function [hold, slow, g] = actuator_gate_state(u, t, faults, g)
% ACTUATOR_GATE_STATE  Memory for actuator_gate: command at onset, slow mode.
%
%   u      the command (sampled at g.dt)
%   hold   per fault-table row: the command when the fault turned on
%   slow   per row: the command rate-limited to the row's magnitude per s
%   g      gate state (gate_state)
%
%#codegen

F = fault_codes();
n = size(faults, 1);
for r = 1:n
    tg = faults(r, 1);
    if tg < F.T_GIMBAL_P || tg > F.T_RCS
        continue
    end
    switch tg
        case F.T_THROTTLE, k = 1;
        case F.T_GIMBAL_P, k = 2;
        case F.T_GIMBAL_Y, k = 3;
        otherwise,         k = 4;
    end
    if t < faults(r, 3) || t >= faults(r, 4)
        g.held_u(r) = false;
        g.u_value(r) = u(k);
        g.u_last(r) = u(k);
        continue
    end
    if ~g.held_u(r)
        g.held_u(r) = true;          % u_value keeps the last pre-onset command
    end
    step = abs(faults(r, 5)) * g.dt;
    g.u_last(r) = g.u_last(r) + min(max(u(k) - g.u_last(r), -step), step);
end
hold = g.u_value;
slow = g.u_last;
end

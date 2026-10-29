function u = actuator_gate(u, t, faults, max_gimbal, hold, slow)
% ACTUATOR_GATE  Apply the fault table's actuator faults to the command.
%
%   u           [thrust; gimbal pitch; gimbal yaw; rcs] going to the
%               actuators (after FDIR's command select)
%   max_gimbal  gimbal limit (rad), for hardover
%   hold, slow  per fault-table row: the command at onset and the
%               rate-limited command (actuator_gate_state)
%
% Stateless, so it runs at the controller's continuous rate: with no fault
% active it returns u untouched and the flight is exactly the nominal one
% (a sampled gate changed the thrust path's If/Merge cutoff logic). The
% memory a fault needs lives in actuator_gate_state, sampled separately.
% The gate sits between the flight computer and the actuators, so the
% command FDIR sent and the actuator feedback differ, which is what the
% gimbal-tracking and RCS-health monitors look for.
%   stuck    hold the command at onset (throttle: magnitude N if > 0)
%   slow     command rate-limited to magnitude per second
%   hardover gimbal pinned at sign(magnitude) * max_gimbal
%   bias     + magnitude
%   dead     zero (throttle: no thrust; rcs: no torque)
%   reduced  scaled by magnitude
%
%#codegen

F = fault_codes();
for r = 1:size(faults, 1)
    tg = faults(r, 1); md = faults(r, 2);
    if tg < F.T_GIMBAL_P || tg > F.T_RCS || t < faults(r, 3) || t >= faults(r, 4)
        continue
    end
    k = actuatorIndex(tg, F);
    mag = faults(r, 5);
    switch md
        case F.M_STUCK
            if tg == F.T_THROTTLE && mag > 0
                u(k) = mag;
            else
                u(k) = hold(r);
            end
        case F.M_SLOW
            u(k) = slow(r);
        case F.M_HARDOVER
            u(k) = sign(mag) * max_gimbal;
        case F.M_BIAS
            u(k) = u(k) + mag;
        case F.M_DEAD
            u(k) = 0;
        case F.M_REDUCED
            u(k) = mag * u(k);
        otherwise
    end
end
end

function k = actuatorIndex(tg, F)
% u index: thrust 1, gimbal pitch 2, gimbal yaw 3, rcs 4
switch tg
    case F.T_THROTTLE, k = 1;
    case F.T_GIMBAL_P, k = 2;
    case F.T_GIMBAL_Y, k = 3;
    otherwise,         k = 4;
end
end

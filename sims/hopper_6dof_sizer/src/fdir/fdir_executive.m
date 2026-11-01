function [y, e] = fdir_executive(raw, metric, u, cat, h_cut, e)
% FDIR_EXECUTIVE  One frame of the FDIR executive (CoDR 15-18, PDR 78-79).
%
%   raw     N x 1 logical, each monitor's pass/fail this frame (catalog order)
%   metric  N x 1 double, the monitor statistic for CUSUM rows
%   u       .t  .phase (fdir_modes)  .alt (m)  .tw (thrust / weight)
%           .reset (ground reset: clears faults on the pad or after landing)
%   cat     fault table (fdir_load_catalog)
%   h_cut   safe height for an in-flight R4 (m)
%   e       state (fdir_exec_state), returned updated
%
%   y.status     N x 1 uint8: 0 INACTIVE, 1 PENDING, 2 LATCHED
%   y.level      active response level: -1 none, 0..4 = R0..R4. The maximum
%                over latched faults; it only escalates (arbitration).
%   y.degrade    an R1 fault is latched (R0 and R1 stack with any mode)
%   y.mask       R0 target bitmask (fdir_targets), sources to exclude
%   y.no_arm     an IBIT fault is latched: hold on the pad
%   y.uncontrollable  attitude envelope latched above h_cut
%   y.first_fault, y.first_fault_t   first fault to latch and when
%
% Per fault: monitors -> persistence (N-of-M, consecutive or CUSUM) ->
% isolation (the fault's target) -> latched flag -> response table. A
% latched fault escalates to escalate_to after escalate_after_s. Faults are
% evaluated only in the phases of their mask and never when inhibited.
%
%#codegen

K = fdir_catalog_cols();
n = size(cat, 1);
H = size(e.history, 2);
INACTIVE = uint8(0); PENDING = uint8(1); LATCHED = uint8(2);

on_ground = u.phase == uint8(1) || u.phase == uint8(6);   % PAD or LANDED
if u.reset && on_ground
    keep = e;
    e = fdir_exec_state(n);
    e.fail_count = keep.fail_count;      % flight totals survive a reset (logged)
end

level = int8(-1);
degrade = false;
mask = 0;
no_arm = false;
uncontrollable = false;

for i = 1:n
    if cat(i, K.INHIBITED) > 0
        continue
    end
    active_phase = bitand(uint8(cat(i, K.PHASE_MASK)), bitshift(uint8(1), double(u.phase) - 1)) > 0;
    fail = logical(raw(i));

    if active_phase && e.status(i) ~= LATCHED
        e.history(i, :) = [e.history(i, 2:end), fail];
        if fail
            e.consec(i) = e.consec(i) + 1;
            e.fail_count(i) = e.fail_count(i) + 1;
        else
            e.consec(i) = 0;
        end
        e.cusum(i) = max(0, e.cusum(i) + metric(i));

        mwin = min(cat(i, K.M), H);
        recent = 0;
        for j = 1:H
            if j > H - mwin && e.history(i, j), recent = recent + 1; end
        end
        switch cat(i, K.FILTER)
            case K.F_N_OF_M,       trip = recent >= cat(i, K.N);
            case K.F_CONSECUTIVE,  trip = e.consec(i) >= cat(i, K.N);
            otherwise,             trip = e.cusum(i) > cat(i, K.CUSUM);
        end
        if trip
            e.status(i) = LATCHED;
            e.first_t(i) = u.t;
            if e.first_fault == 0
                e.first_fault = i;
                e.first_fault_t = u.t;
            end
        elseif fail || recent > 0
            e.status(i) = PENDING;
        else
            e.status(i) = INACTIVE;
        end
    elseif ~active_phase && e.status(i) ~= LATCHED
        e.history(i, :) = false;
        e.consec(i) = 0;
        e.status(i) = INACTIVE;
    elseif e.status(i) == LATCHED && cat(i, K.CLEAR) == K.C_AUTO_CLEAR && active_phase
        % auto-clearing faults release after a full window of passes
        e.history(i, :) = [e.history(i, 2:end), fail];
        mwin = min(cat(i, K.M), H);
        recent = 0;
        for j = 1:H
            if j > H - mwin && e.history(i, j), recent = recent + 1; end
        end
        if recent == 0
            e.status(i) = INACTIVE;
            e.consec(i) = 0;
            e.cusum(i) = 0;
        end
    end

    if e.status(i) ~= LATCHED
        continue
    end

    % ---- response table ----
    r = cat(i, K.RESPONSE);
    if r == K.R_NO_ARM
        no_arm = true;
        continue
    elseif r == K.R_R4_ABOVE_HCUT
        if u.alt > h_cut
            r = 4; uncontrollable = true;
        else
            r = 2;
        end
    elseif r == K.R_STUCK_THRUST
        if u.tw > 1.1, r = 4; else, r = 2; end
    end
    if cat(i, K.ESC_TO) >= 0 && u.t - e.first_t(i) >= cat(i, K.ESC_AFTER)
        r = max(r, cat(i, K.ESC_TO));
    end
    if r == 1, degrade = true; end
    if r == 0
        mask = bitor(mask, bitshift(1, cat(i, K.TARGET) - 1));
    end
    level = max(level, int8(r));
end

e.level = max(e.level, level);      % the ladder only escalates
y.status = e.status;
y.level = e.level;
y.degrade = degrade;
y.mask = mask;
y.no_arm = no_arm;
y.uncontrollable = uncontrollable;
y.first_fault = e.first_fault;
y.first_fault_t = e.first_fault_t;
end

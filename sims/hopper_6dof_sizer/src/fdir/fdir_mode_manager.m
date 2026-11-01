function [out, s] = fdir_mode_manager(u, p, s)
% FDIR_MODE_MANAGER  One frame of the GNC PDR flight-mode state machine.
%
%   u  this frame's inputs
%        .t                 sim time (s)
%        .alt               height above the pad (m, up +)
%        .vz                vertical velocity (m/s, up +)
%        .thrust_ratio      delivered thrust / weight
%        .response          active FDIR response level (fdir_modes R0..R4)
%        .uncontrollable    loss-of-control flag (R4 allowed at any height)
%        .touchdown         touchdown confirmed (fdir_touchdown)
%        .thrust_cut        nominal landing thrust cut commanded
%        .no_arm            an IBIT check failed: hold before ARMED
%        .cmd               ground commands, logical: init, arm, start,
%                           launch, abort, terminate, reset
%   p  parameters (fdir_mode_params)
%   s  state (fdir_mode_state), returned updated
%
%   out.mode       mode in effect this frame (a transition decided now
%                  takes effect next frame, PDR slide 78)
%   out.phase      that mode's phase (fault-table phase mask)
%   out.abort      R2/R3 abort-to-land guidance active
%   out.nav_abort  R3 (runs as R2 until NAVIGATION provides the shadow
%                  estimator; flagged so the run is marked)
%   out.safed      R4 cutoff and safing active
%   out.r4_permitted  R4 would be allowed now (pad, below h_cut, or
%                  uncontrollable)
%
% Contingencies (PDR slides 79, 84, 86): abort (R2) and terminate (R4)
% ground commands are always accepted. An automatic R4 in flight waits,
% flying the abort, until the vehicle is below the safe height h_cut or
% uncontrollable; on the pad it is always allowed. R2/R3 exit only on
% confirmed touchdown. R0/R1 never change the mode.
%
%#codegen

M = fdir_modes();
mode = s.mode;
out.mode = mode;
out.phase = M.MODE_PHASE(mode);

in_pad_modes = mode <= M.PAD_CHECKOUT;
airborne_modes = mode >= M.ASCENT && mode <= M.TERMINAL_DESCENT;
r4_permitted = in_pad_modes || u.alt <= p.h_cut || u.uncontrollable;
out.r4_permitted = r4_permitted;

next = mode;

% ---- Contingencies first: they override the nominal sequence ----
r4_requested = u.response >= M.R4 || u.cmd.terminate;
if mode ~= M.SAFED && mode ~= M.LANDED && r4_requested && (r4_permitted || u.cmd.terminate)
    next = M.SAFED;
    s.r4_pending = false;
elseif mode == M.SAFED || mode == M.LANDED
    % stay until a ground reset
elseif mode >= M.ENGINE_START && in_pad_modes && (u.response >= M.R2 || u.cmd.abort)
    next = M.SAFED;                      % on the pad an abort is a shutdown
elseif (airborne_modes || mode == M.ABORT_LAND || mode == M.NAV_ABORT) && ...
        (u.response >= M.R2 || u.cmd.abort || r4_requested)
    s.r4_pending = r4_requested;   % R4 not yet permitted: fly the abort meanwhile
    if u.response == M.R3 && ~r4_requested
        next = M.NAV_ABORT;
    elseif mode ~= M.NAV_ABORT
        next = M.ABORT_LAND;
    end
end

% ---- Nominal sequence (only when no contingency moved the mode) ----
if next == mode
    switch mode
        case M.SAFE
            if u.cmd.init || p.auto_sequence, next = M.INITIALIZE; end
        case M.INITIALIZE
            next = M.READY;                      % initialization is instantaneous in the sim
        case M.READY
            if ~u.no_arm && (u.cmd.arm || p.auto_sequence), next = M.ARMED; end
        case M.ARMED
            if u.cmd.start || p.auto_sequence, next = M.ENGINE_START; end
        case M.ENGINE_START
            if u.thrust_ratio >= p.engine_stable_ratio, next = M.PAD_CHECKOUT; end
        case M.PAD_CHECKOUT
            if (u.cmd.launch || p.auto_sequence) && u.alt > p.liftoff_alt && u.vz > 0
                next = M.ASCENT;
            end
        case M.ASCENT
            if abs(u.alt - p.hover_alt) <= p.hover_capture_alt && abs(u.vz) <= p.hover_capture_vel
                next = M.HOVER;
                s.hover_start = u.t;
            end
        case M.HOVER
            if u.t - s.hover_start >= p.hover_time
                next = M.DESCENT;
            end
        case M.DESCENT
            if u.alt <= p.terminal_alt, next = M.TERMINAL_DESCENT; end
        case M.TERMINAL_DESCENT
            if u.thrust_cut, next = M.TOUCHDOWN_CONFIRM; end
        case M.TOUCHDOWN_CONFIRM
            if u.touchdown, next = M.LANDED; end
        case {M.ABORT_LAND, M.NAV_ABORT}
            if u.touchdown, next = M.LANDED; end
        otherwise
    end
end

% ---- Ground reset after landing or safing ----
if (mode == M.LANDED || mode == M.SAFED) && u.cmd.reset
    next = M.SAFE;
    s.hover_start = -1.0;
    s.r4_pending = false;
end

if next ~= mode
    s.mode = next;
    s.mode_entered = u.t;
end

out.abort = mode == M.ABORT_LAND || mode == M.NAV_ABORT;
out.nav_abort = mode == M.NAV_ABORT;
out.safed = mode == M.SAFED;
end

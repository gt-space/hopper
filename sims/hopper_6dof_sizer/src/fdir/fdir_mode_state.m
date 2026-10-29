function s = fdir_mode_state()
% FDIR_MODE_STATE  Initial state for fdir_mode_manager and fdir_touchdown.
%
%#codegen

M = fdir_modes();
s.mode         = M.SAFE;
s.mode_entered = 0.0;     % time the current mode was entered (s)
s.hover_start  = -1.0;    % time HOVER was entered (s), -1 = not yet
s.r4_pending   = false;   % R4 latched but not yet permitted (flying R2 meanwhile)
s.first_mode_change = 0.0;

% touchdown persistence: last M evaluations of "two independent signals"
s.td_history = false(1, 10);
end

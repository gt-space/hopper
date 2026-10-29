function p = fdir_mode_params(IN)
% FDIR_MODE_PARAMS  Mode-manager and touchdown parameters for a vehicle.
%
%   p = fdir_mode_params(IN)   IN from mission_inputs (hopper or CTV)
%
% The hover and terminal-descent limits follow the reference trajectory in
% IN.mission.trajectory, so they scale with the vehicle. IN.fdir.modes
% overrides any field. Values not set by a requirement are provisional:
% the PDR sets them from Monte Carlo statistics (slide 82, "Open").

T = IN.mission.trajectory;

p.auto_sequence       = true;    % sim authorizes init / arm / start / launch itself
p.engine_stable_ratio = 0.8;     % ENGINE START -> PAD CHECKOUT: thrust / weight
p.liftoff_alt         = 0.02;    % m, PAD CHECKOUT -> ASCENT once above this and climbing
p.hover_alt           = T.apex;  % m
p.hover_capture_alt   = max(0.5, 0.02 * T.apex);   % m
p.hover_capture_vel   = 0.3;     % m/s
p.hover_time          = T.t_hover;                 % s
p.terminal_alt        = 2.5 * T.z_td;              % m, DESCENT -> TERMINAL DESCENT
p.h_cut               = 1.5;     % m, safe height for an in-flight R4 (PDR: ~1-2 m, set by legs)

% touchdown (PDR slide 85)
p.contact_alt    = 0.05;   % m, lidar median at or below this
p.impulse_thresh = 3.0;    % m/s^2, specific-force jump at leg contact
p.vz_still       = 0.1;    % m/s
p.td_n           = 3;      % confirm when 3 of the last 5 frames see two signals
p.td_m           = 5;

if isfield(IN, 'fdir') && isfield(IN.fdir, 'modes')
    for f = fieldnames(IN.fdir.modes)'
        p.(f{1}) = IN.fdir.modes.(f{1});
    end
end
end

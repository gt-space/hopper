function [confirmed, n_signals, s] = fdir_touchdown(u, p, mode, s)
% FDIR_TOUCHDOWN  Touchdown confirmation (GNC PDR slide 85).
%
%   u.lidar_median  median lidar range (m)
%   u.accel_impulse magnitude of the specific-force jump this frame (m/s^2)
%   u.vz            vertical velocity (m/s, up +)
%   u.thrust_cut    thrust cut commanded
%
% Touchdown needs two independent signals, held N of the last M frames.
% Detection is armed only in TERMINAL DESCENT, TOUCHDOWN CONFIRM or an
% abort, so nothing at altitude can trigger it. One signal alone does
% nothing.
%
%#codegen

M = fdir_modes();
armed = mode == M.TERMINAL_DESCENT || mode == M.TOUCHDOWN_CONFIRM || ...
        mode == M.ABORT_LAND || mode == M.NAV_ABORT;

sig = [u.lidar_median <= p.contact_alt, ...
       u.accel_impulse >= p.impulse_thresh, ...
       abs(u.vz) <= p.vz_still, ...
       logical(u.thrust_cut)];
n_signals = uint8(sum(sig));

hit = armed && n_signals >= 2;
s.td_history = [s.td_history(2:end), hit];
% N of the last M frames, as a fixed-bound loop so it compiles in a MATLAB
% Function block without dynamic sizing
m_max = numel(s.td_history);
count = 0;
for j = 1:m_max
    if j > m_max - p.td_m && s.td_history(j)
        count = count + 1;
    end
end
confirmed = armed && count >= p.td_n;
end

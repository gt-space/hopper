function m = baffleSimMetrics(so, t_ref, alt_ref, alt_split)
% Altitude / pitch metrics for the baffle sim sweep (PostSimFcn).
%   pitch_err       - max |pitch - 90| over the flight, deg
%   pitch_err_high  - same, while altitude >= alt_split (default 1 m)
%   pitch_err_low   - same, on the final descent below alt_split
%                     (landing transient)
%   alt_err         - max |altitude - reference altitude| where the
%                     reference is defined, m
%   max_alt         - peak altitude, m

if nargin < 4, alt_split = 1; end

m.ok = isempty(so.ErrorMessage);
m.error = so.ErrorMessage;
m.pitch_err = NaN; m.pitch_err_high = NaN; m.pitch_err_low = NaN;
m.alt_err = NaN; m.max_alt = NaN; m.t_end = NaN;
if ~m.ok
    return
end

t     = so.z.Time;
alt   = -so.z.Data;
pitch = so.pitch.Data;
perr  = abs(pitch - 90);

alt_r = interp1(t_ref(:), alt_ref(:), t, 'linear', NaN);

% final descent below alt_split (after apogee)
[~, i_apo] = max(alt);
i_low = find((1:numel(t))' > i_apo & alt < alt_split, 1);
if isempty(i_low), i_low = numel(t) + 1; end

m.pitch_err      = max(perr);
m.pitch_err_high = max(perr(1:i_low-1));
m.pitch_err_low  = max([perr(i_low:end); 0]);
m.alt_err        = max(abs(alt - alt_r), [], 'omitnan');
m.max_alt        = max(alt);
m.t_end          = t(end);
end

function S = tvcEmInit() %#codegen
%TVCEMINIT  Zero state for tvcEmUpdate (two axes: pitch, yaw).
S.i  = zeros(2, 1);      % A, winding current
S.x  = zeros(2, 1);      % m, rod position (positive toward +delta)
S.xd = zeros(2, 1);      % m/s, rod speed
S.I  = zeros(2, 1);      % rad s, servo integrator
S.dg = zeros(2, 1);      % rad, gimbal angle (after backlash)
S.V  = zeros(2, 1);      % V, held servo output
S.buf = zeros(101, 2);   % command delay line, one row per servo period
end

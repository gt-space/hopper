function S = tvcAxisSim(P, ~, t_end, cmdFcn, opts)
%TVCAXISSIM  Bench-style run of one TVC axis on the electromechanical model.
%   S = tvcAxisSim(P, [], 1.0, @(t) deg2rad(2)*(t>=0))
%
%   Steps tvcEmUpdate (the same code the Simulink TVC_Actuator variant
%   runs) once per servo period with the pitch axis commanded and yaw held
%   at zero, at constant thrust P.thrust. Outputs are sampled at the servo
%   period (1 ms).
%
%   opts.e        thrust-line offset loading the axis (m), default P.e_nom
%   opts.coulomb, opts.backlash  logical, default true

arguments
    P struct
    ~
    t_end (1,1) double
    cmdFcn function_handle
    opts.e (1,1) double = NaN
    opts.coulomb (1,1) logical = true
    opts.backlash (1,1) logical = true
    opts.decimate (1,1) double = 1   % unused, kept for older callers
    opts.accel_lim (1,1) double = Inf  % rad/s^2, drive acceleration ramp; Inf = off
    opts.rate_lim  (1,1) double = Inf  % rad/s,   drive speed clamp; Inf = off
end
e = opts.e;  if isnan(e), e = P.e_nom; end
T = tvcEmParams(P, Offset=[e; e], Coulomb=opts.coulomb, Backlash=opts.backlash, MaxAngle=deg2rad(20), AccelLimit=opts.accel_lim, RateLimit=opts.rate_lim);
N = round(t_end / T.dt_c);
St = tvcEmInit();
S.t = (0:N)' * T.dt_c;
S.delta = zeros(N + 1, 1);  S.cmd = S.delta;  S.rate = S.delta;  S.i = S.delta;  S.V = S.delta;
for n = 1:N
    c = cmdFcn(S.t(n));
    [St, d, r, i] = tvcEmUpdate(St, [c; 0], P.thrust, T);
    S.cmd(n) = c;  S.delta(n + 1) = d(1);  S.rate(n + 1) = r(1);  S.i(n + 1) = i(1);  S.V(n + 1) = St.V(1);
end
S.cmd(end) = cmdFcn(S.t(end));
S.F = P.KF * S.i;
S.accel = gradient(S.rate, T.dt_c);
end

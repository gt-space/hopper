function T = tvcEmParams(P, opts)
%TVCEMPARAMS  Numeric parameter struct for the electromechanical TVC model.
%   T = tvcEmParams(tvcActuatorParams('AG9113', Layout='recommended'))
%   T = tvcEmParams(P, Latency=0.01, Offset=[0.01; 0.01], MaxAngle=deg2rad(8))
%
%   T is what tvcEmInit / tvcEmUpdate read: the Simulink TVC_Actuator
%   variant (Actuators_TVC_RCS, TVC_MODEL = ag9113 / ag9313) gets it as
%   the workspace variable TVCEM, and tvcAxisSim uses it for the
%   standalone trade, so both run the same code. Numeric fields only
%   (code generation).
%
%   Offset is the thrust-line offset (m) whose moment T*e loads the
%   [pitch; yaw] actuators, opposing +motion (conservative constant).

arguments
    P struct
    opts.Latency (1,1) double = P.lat
    opts.Offset (2,1) double = [P.e_nom; P.e_nom]
    opts.MaxAngle (1,1) double = deg2rad(8)
    opts.Coulomb (1,1) logical = true
    opts.Backlash (1,1) logical = true
    opts.AccelLimit (1,1) double = Inf   % rad/s^2, drive acceleration ramp; Inf = off
    opts.RateLimit  (1,1) double = Inf   % rad/s,   drive speed limit; Inf = off
end
K = tvcKinematics(P, -1:1, -1:1);   % only K.axis is used
d = K.axis.delta;  x = K.axis.dL;  r = K.axis.r;
s = sign(interp1(d, r, 0));          % rod coordinate positive toward +delta
x = s * x;  r = s * r;
m = abs(d) <= deg2rad(16);

T.Vbus = P.Vbus;  T.Imax = P.Imax;  T.Ke = P.Ke;  T.KF = P.KF;  T.R0 = P.R0;
T.L = P.L;  T.m_r = P.m_r;
T.Fc = P.Fc * opts.Coulomb;  T.backlash = P.backlash * opts.Backlash;
T.Jg = P.Jg;  T.kh = P.kh;  T.tauf = P.tauf;
T.Kp = P.Kp;  T.Ki = P.Ki;  T.Kd = P.Kd;
T.e = opts.Offset;
T.pD = polyfit(x(m), d(m), 5);       % delta(x)
T.pR = polyfit(d(m), r(m), 5);       % r(delta) = dx/ddelta
T.dt_c = P.dt_c;                     % s, built-in servo period (= model step)
T.n_sub = 100;                       % integration sub-steps per servo period
T.lat_steps = round(opts.Latency / P.dt_c);
T.max_angle = opts.MaxAngle;         % rad, commanded-angle limit
T.accel_lim = opts.AccelLimit;       % rad/s^2, setpoint acceleration ramp (Inf = off)
T.rate_lim  = opts.RateLimit;        % rad/s,   setpoint speed clamp (Inf = off)
if T.lat_steps > 100
    error('tvcEmParams:latency', 'Latency above 100 servo periods is not supported.');
end
end

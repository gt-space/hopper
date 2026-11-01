out = sim('throttle_valve_simulink');

Pout = out.Pout_log.signals.values;
t    = out.Pout_log.time;

Pout = Pout(:);
t    = t(:);

info_outer = stepinfo(Pout, t);
%                     ^signal       ^time     ^final expected value (your P_set = 100)
info_inner = stepinfo(out.theta_log, out.tout);
%                     ^signal        ^time     ^final expected value (max theta = 90°)

BW_outer = 0.35 / info_outer.RiseTime;
BW_inner = 0.35 / info_inner.RiseTime;
fprintf('Outer BW: %.4f rad/s\n', BW_outer);
fprintf('Inner BW: %.4f rad/s\n', BW_inner);
fprintf('Ratio (inner/outer): %.11fx -- target is 5-10x\n', BW_inner / BW_outer);
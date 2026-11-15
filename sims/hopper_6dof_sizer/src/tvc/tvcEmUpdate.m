function [S, delta, rate, cur] = tvcEmUpdate(S, cmd, thrust, T) %#codegen
%TVCEMUPDATE  One servo period of the electromechanical TVC model, 2 axes.
%   [S, delta] = tvcEmUpdate(S, cmd, thrust, T)
%     cmd     [pitch; yaw] commanded gimbal angle (rad)
%     thrust  engine thrust (N), loads the actuators through T.e
%     T       tvcEmParams struct;  S from tvcEmInit
%     delta   gimbal angles at the end of the period (rad)
%     rate    gimbal rates (rad/s);  cur  winding currents (A)
%
%   GNC PDR slide 170 (simplified proposed model) plus the TBD terms:
%     V      = sat(Kp e + Ki int e - Kd deltadot, +/-Vbus), held for T.dt_c,
%              e = cmd(t - lat) - delta_motor, integrator frozen while saturated
%     L di/dt = V - R0 i - Ke xdot,  |i| <= Imax
%     m_eq xddot = KF i - Fc sgn xdot - tau_h / r(delta),  m_eq = m_r + Jg/r^2
%     tau_h  = kh delta + tauf sgn deltadot + thrust * e
%     delta  = IK^-1(x) from the sizer kinematics (T.pD, T.pR), with a
%              backlash band between the motor side and the gimbal
%   Integrated with T.n_sub semi-implicit Euler steps per servo period.

dt = T.dt_c / T.n_sub;
delta = zeros(2, 1);  rate = zeros(2, 1);  cur = zeros(2, 1);
for a = 1:2
    % ---- command delay line and limit ----
    S.buf(2:end, a) = S.buf(1:end-1, a);
    S.buf(1, a) = min(max(cmd(a), -T.max_angle), T.max_angle);
    c = S.buf(T.lat_steps + 1, a);

    % ---- built-in servo (once per period) ----
    dm = polyval(T.pD, S.x(a));
    r  = polyval(T.pR, dm);
    err = c - dm;
    Vu = T.Kp * err + T.Ki * S.I(a) - T.Kd * S.xd(a) / r;
    V = min(max(Vu, -T.Vbus), T.Vbus);
    if V == Vu
        S.I(a) = S.I(a) + err * T.dt_c;
    end
    S.V(a) = V;

    % ---- motor, screw and gimbal ----
    tauT = thrust * T.e(a);
    i = S.i(a);  x = S.x(a);  xd = S.xd(a);  dg = S.dg(a);
    for n = 1:T.n_sub
        dm = polyval(T.pD, x);
        r  = polyval(T.pR, dm);
        if T.backlash > 0
            if dm - dg > T.backlash / 2
                dg = dm - T.backlash / 2;
            elseif dg - dm > T.backlash / 2
                dg = dm + T.backlash / 2;
            end
        else
            dg = dm;
        end
        tau_h = T.kh * dg + T.tauf * tanh(xd / r / 1e-3) + tauT;
        m_eq = T.m_r + T.Jg / r^2;
        i  = i + dt / T.L * (V - T.R0 * i - T.Ke * xd);
        i  = min(max(i, -T.Imax), T.Imax);
        xd = xd + dt / m_eq * (T.KF * i - T.Fc * tanh(xd / 1e-4) - tau_h / r);
        x  = x + dt * xd;
    end
    S.i(a) = i;  S.x(a) = x;  S.xd(a) = xd;  S.dg(a) = dg;
    delta(a) = dg;
    rate(a) = xd / polyval(T.pR, polyval(T.pD, x));
    cur(a) = i;
end
end

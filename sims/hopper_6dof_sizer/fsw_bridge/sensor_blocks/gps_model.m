function gps = gps_model(p_n, v_n, C_bn, w_b, thrust, SENS)
%#codegen
% GPS measurement model (Notion: GPS Measurement Model), u-blox ZED-F9P.
%   p_n, v_n  CG position (m) and velocity (m/s), NED from the pad
%   C_bn      NED -> body rotation
%   w_b       body angular rate (rad/s)
%   thrust    delivered thrust (N); engine on raises the dropout chance
%   gps       [lat deg; lon deg; alt m; vN; vE; vD m/s; has_fix; num_sats]
%             (luna GpsState fields)
%
%   p = p_cg + Cbn' r_b + b_slow + b_fast + eta_p
%   v = v_cg + Cbn' (w x r_b) + eta_v
%   Drifts are unit Gauss-Markov states scaled by the mode's sigmas. Each fix
%   is converted to lat/lon/height, quantized, and delivered `latency` later.
%   Runs every `tick` seconds so fix times and latency land on the grid.
%   Not modelled yet: pad multipath, RTK fixed/float transitions.
P = SENS.gps;
m = P.mode + 1;
persistent s zs zf out pend pend_at k fix
if isempty(s)
    s = rng_seed(SENS.seed, 4);
    [zs, s] = rng_gauss(s, 3);
    [zf, s] = rng_gauss(s, 3);
    out = zeros(8, 1);          % nothing received yet: has_fix = 0
    pend = zeros(8, 1);
    pend_at = -1;
    k = 0;
    fix = true;
end

Tg = 1 / P.rate;
n_fix = round(Tg / P.tick);
n_lat = round(P.latency / P.tick);

if mod(k, n_fix) == 0
    % Drift states advance once per fix
    phs = exp(-Tg / P.tau_slow(m));
    phf = exp(-Tg / P.tau_fast(m));
    [n, s] = rng_gauss(s, 3); zs = phs * zs + sqrt(1 - phs^2) * n;
    [n, s] = rng_gauss(s, 3); zf = phf * zf + sqrt(1 - phf^2) * n;

    % Fix lost / regained
    [u, s] = rng_uniform(s);
    if fix
        p_loss = P.p_loss(1);
        if thrust > 0
            p_loss = P.p_loss(2);
        end
        fix = u >= p_loss;
    else
        fix = u < 1 - exp(-Tg / P.T_reacq);
    end

    [np, s] = rng_gauss(s, 3);
    [nv, s] = rng_gauss(s, 3);
    if fix
        dH = P.dop_scale(1);
        dV = P.dop_scale(2);
        C_nb = C_bn';
        pos = p_n(:) + C_nb * P.r_b ...
            + [P.sH_slow(m) * dH; P.sH_slow(m) * dH; P.sV_slow(m) * dV] .* zs ...
            + [P.sH_fast(m) * dH; P.sH_fast(m) * dH; P.sV_fast(m) * dV] .* zf ...
            + [P.sH_white(m) * dH; P.sH_white(m) * dH; P.sV_white(m) * dV] .* np;
        vel = v_n(:) + C_nb * cross(w_b(:), P.r_b) + [P.svH * dH; P.svH * dH; P.svD * dV] .* nv;

        % NED metres from the pad -> geodetic (flat Earth over the flight)
        a = 6378137;
        e2 = 0.00669438;
        lat0 = P.lat0 * pi / 180;
        sl2 = sin(lat0)^2;
        N0 = a / sqrt(1 - e2 * sl2);
        RM = a * (1 - e2) / (1 - e2 * sl2)^1.5;
        lat = P.lat0 + pos(1) / (RM + P.h0) * 180 / pi;
        lon = P.lon0 + pos(2) / ((N0 + P.h0) * cos(lat0)) * 180 / pi;
        alt = P.h0 - pos(3);

        q = P.dq_deg(m);
        qa = P.dq_alt(m);
        pend = [q * round(lat / q); q * round(lon / q); qa * round(alt / qa); ...
                P.dq_vel * round(vel / P.dq_vel); 1; P.num_sats];
    else
        pend = [out(1:6); 0; 0];  % no new solution; last values held
    end
    pend_at = k + n_lat;
end

if k == pend_at
    out = pend;
end
k = k + 1;
gps = out;
end

% --- deterministic RNG shared by the sensor blocks (xorshift32 + Box-Muller)
function s = rng_seed(seed, stream)
s = uint32(mod(double(seed) * 2654435761 + double(stream) * 97531 + 12345, 4294967296));
if s == 0
    s = uint32(1);
end
end

function [u, s] = rng_uniform(s)
s = bitxor(s, bitshift(s, 13));
s = bitxor(s, bitshift(s, -17));
s = bitxor(s, bitshift(s, 5));
u = (double(s) + 0.5) / 4294967296;
end

function [n, s] = rng_gauss(s, k)
n = zeros(k, 1);
for i = 1:k
    [u1, s] = rng_uniform(s);
    [u2, s] = rng_uniform(s);
    n(i) = sqrt(-2 * log(u1)) * cos(2 * pi * u2);
end
end

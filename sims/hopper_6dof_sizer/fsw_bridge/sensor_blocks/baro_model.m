function baro = baro_model(p_n, v_n, C_bn, thrust, SENS)
%#codegen
% Barometer measurement model (Notion: Barometer Model), MS5611.
%   p_n    CG position NED (m), ground at D = 0
%   v_n    CG velocity NED (m/s)
%   C_bn   NED -> body rotation
%   thrust delivered thrust (N); > 0 means the engine is on
%   baro   [pressure Pa; temperature degC]
%
%   p = p_atm(h) + dp_ground_effect(h) + dp_dynamic + beta + eta, quantized
%   Tare/zeroing is left to the flight software.
P = SENS.baro;
persistent s b
if isempty(s)
    s = rng_seed(SENS.seed, 3);
    [n, s] = rng_gauss(s, 1); b = P.b0_sd * n;
end

% Height of the sensor port above the pad
r_n = C_bn' * P.r_b;
h = -(p_n(3) + r_n(3));

% Standard atmosphere from the pad conditions
T = P.T_pad - P.L * h;
p_atm = P.p_pad * (T / P.T_pad)^(P.g0 / (P.L * P.R));

% Exhaust ground effect, only near the pad with the engine on
dp_ge = 0;
if thrust > 0 && h < P.ge_D
    dp_ge = P.ge_sign * P.ge_A * (P.ge_D - max(h, 0))^2;
end

% Part of the dynamic pressure reaching the port
rho = p_atm / (P.R * T);
v_rel = v_n(:) - P.wind_n;
dp_dyn = P.k_port * 0.5 * rho * (v_rel' * v_rel);

% Gauss-Markov bias drift and white noise
phi = exp(-P.Ts / P.b_tau);
[n, s] = rng_gauss(s, 2);
b = phi * b + sqrt(1 - phi^2) * P.b_sd * n(1);
p = p_atm + dp_ge + dp_dyn + b + P.noise_sd * n(2);

baro = [P.dq * round(p / P.dq); T - 273.15];
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

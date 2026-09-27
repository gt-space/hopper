function imu = imu_model(F_b, W_b, mass, w_b, SENS)
%#codegen
% IMU measurement model (Notion: IMU Measurement Model), ADIS16500.
%   F_b  total force on the vehicle, body axes (N)
%   W_b  weight, body axes (N)
%   mass vehicle mass (kg)
%   w_b  body angular rate (rad/s)
%   imu  [accel (3) m/s^2; gyro (3) deg/s], sensor axes
%
%   f = (I+Ma)(I+Sa) Rbs f_true + ba + na
%   w = (I+Mg)(I+Sg) Rbs w_true + Gg Rbs f_true + bg + ng
%   Biases = turn-on offset (drawn once) + first-order Gauss-Markov drift.
P = SENS.imu;
persistent s ba0 bg0 bad bgd Sa Sg Ma Mg
if isempty(s)
    s = rng_seed(SENS.seed, 1);
    [n, s] = rng_gauss(s, 3); ba0 = P.acc_b0_sd * n;
    [n, s] = rng_gauss(s, 3); bg0 = P.gyro_b0_sd * n;
    [n, s] = rng_gauss(s, 3); bad = P.acc_b_sd * n;
    [n, s] = rng_gauss(s, 3); bgd = P.gyro_b_sd * n;
    [n, s] = rng_gauss(s, 3); Sa = diag(1 + P.sf_sd * n);
    [n, s] = rng_gauss(s, 3); Sg = diag(1 + P.sf_sd * n);
    [n, s] = rng_gauss(s, 6); Ma = eye(3) + P.misalign_sd * [0 n(1) n(2); n(3) 0 n(4); n(5) n(6) 0];
    [n, s] = rng_gauss(s, 6); Mg = eye(3) + P.misalign_sd * [0 n(1) n(2); n(3) 0 n(4); n(5) n(6) 0];
end

% Specific force at the CG, plus centripetal term at the sensor location
f_true = (F_b(:) - W_b(:)) / mass + cross(w_b(:), cross(w_b(:), P.r_b));
f_s = P.R_bs * f_true;
w_s = P.R_bs * w_b(:);

% Gauss-Markov bias drift
phi_a = exp(-P.Ts / P.acc_b_tau);
phi_g = exp(-P.Ts / P.gyro_b_tau);
[n, s] = rng_gauss(s, 3); bad = phi_a * bad + sqrt(1 - phi_a^2) * P.acc_b_sd * n;
[n, s] = rng_gauss(s, 3); bgd = phi_g * bgd + sqrt(1 - phi_g^2) * P.gyro_b_sd * n;

% White noise, per-sample sd from the noise density
[na, s] = rng_gauss(s, 3);
[ng, s] = rng_gauss(s, 3);
f_meas = Ma * Sa * f_s + ba0 + bad + P.acc_nd / sqrt(P.Ts) * na;
w_meas = Mg * Sg * w_s + P.g_sens * f_s + bg0 + bgd + P.gyro_nd / sqrt(P.Ts) * ng;

imu = [f_meas; w_meas * (180 / pi)];
end

% --- deterministic RNG shared by the sensor blocks (xorshift32 + Box-Muller)
% Same numbers in Simulink and in generated C, seeded per run and per sensor.
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

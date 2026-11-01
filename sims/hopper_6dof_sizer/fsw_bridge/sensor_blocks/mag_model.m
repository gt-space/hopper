function mag = mag_model(C_bn, SENS)
%#codegen
% Magnetometer measurement model (Notion: Magnetometer Model), LIS2MDL.
%   C_bn NED -> body rotation (6DOF DCMbe)
%   mag  magnetic field, sensor axes, Gauss
%
%   z = A Cbn Bn + beta + eta
%   A folds soft iron, scale factor and misalignment; beta is hard iron.
%   Both are drawn once per run (they are what calibration estimates).
P = SENS.mag;
persistent s A beta
if isempty(s)
    s = rng_seed(SENS.seed, 2);
    [n, s] = rng_gauss(s, 9); A = eye(3) + P.A_sd * reshape(n, 3, 3);
    [n, s] = rng_gauss(s, 3); beta = P.beta_sd * n;
end

[n, s] = rng_gauss(s, 3);
mag = A * (C_bn * P.B_n) + beta + P.noise_sd * n;
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

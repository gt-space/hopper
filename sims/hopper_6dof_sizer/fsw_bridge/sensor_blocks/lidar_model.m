function lidar = lidar_model(p_n, C_bn, SENS)
%#codegen
% LiDAR measurement model, 4x Benewake TF-Luna (simplified).
%   p_n    CG position NED (m), ground plane at D = 0
%   C_bn   NED -> body rotation
%   lidar  range along each beam (m), -1 when there is no valid return
%
%   Range from each sensor along its beam to the flat ground, with
%   range-dependent noise, 1 cm steps, the datasheet range window and random
%   dropouts. Stand-in until the ray-traced model on Notion (plume, dust,
%   surface reflectivity) is implemented.
P = SENS.lidar;
persistent s
if isempty(s)
    s = rng_seed(SENS.seed, 5);
end

C_nb = C_bn';
lidar = -ones(4, 1);
for i = 1:4
    pos = p_n(:) + C_nb * P.r_b(:, i);
    d = C_nb * P.d_b(:, i);
    h = -pos(3);
    [u, s] = rng_uniform(s);
    [n, s] = rng_gauss(s, 1);
    if d(3) > 1e-6 && h > 0
        R = h / d(3);
        if R >= P.r_min && R <= P.r_max && u >= P.p_drop
            sd = P.sd_near;
            if R > 3
                sd = P.sd_frac * R;
            end
            lidar(i) = P.dq * round((R + sd * n) / P.dq);
        end
    end
end
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

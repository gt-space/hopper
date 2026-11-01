function [zeta_axial] = axialDampingCalc( ...
    prop_mass_profile, prop_density, nu, tank_r, baffle_w, Nb, tank_h, ...
    baffle_t, slosh_amp)
% Damping ratio of the first axisymmetric (m = 0) slosh mode in a
% cylindrical tank with Nb flat ring baffles.
%
% Energy method, same derivation Miles (1958) used for the lateral mode:
%   zeta = (energy dissipated per cycle) / (4*pi * modal energy)
%
% Mode shape (potential flow), k = 3.8317 / R (1st root of J0' = -J1):
%   eta(r)  = eta_w * J0(k r) / J0(k R)                surface elevation
%   X(r, z) = eta(r) * sinh(k z) / sinh(k h)           vertical displacement
%   E       = 1/2 * rho * g * pi * R^2 * eta_w^2       modal energy
%
% Baffles: quadratic drag on the ring in the vertical oscillating flow,
% Keulegan-Carpenter drag coefficient as used by Miles,
%   C_D = 15 * (U_m T / 2w)^(-1/2),   U_m T = 2*pi*X
%   D   = (4 rho omega^2 / 3) * Int_ring C_D * X^3 dA
% Each point on the ring only dissipates while submerged (Bauer,
% MTP-AERO-62-81); thickness effectiveness from Cole, NASA TN D-3716.
% Applied to the m = 1 mode this reproduces Miles' lateral formula,
% 2.83 e^(-4.6 d/a) alpha^1.5 sqrt(eta_w/a), to within ~9%.
%
% Walls: laminar Stokes boundary layer dissipation on the side wall and
% floor, mean power (rho/2) sqrt(nu omega / 2) |u_slip|^2 per unit area.
%
% The baffle term is independent of g (D and E both scale with omega^2);
% g only enters the wall term through omega.
%
% baffle_t  - baffle thickness [m] (default 0 -> thin plate, no correction)
% slosh_amp - free surface slosh amplitude at the tank wall, eta_w [m]
%             (default 0.1 * tank_r)

if nargin < 8, baffle_t  = 0;            end
if nargin < 9, slosh_amp = 0.1 * tank_r; end

g      = 9.81;
eps_01 = 3.8317;          % 1st axisymmetric mode, root of J1
k      = eps_01 / tank_r;
A_tank = pi * tank_r^2;

% arrays: dim 1 fill level, dim 2 baffle, dim 3 radial point on the ring
h     = max((prop_mass_profile(:) / prop_density) / A_tank, 1e-3);
z_baf = linspace(0, tank_h, Nb);

r_ring   = reshape(linspace(tank_r - baffle_w, tank_r, 50), 1, 1, []);
eta_ring = slosh_amp * abs(besselj(0, k * r_ring) / besselj(0, k * tank_r));

omega   = sqrt(g * k * tanh(k * h));
E_modal = 0.5 * prop_density * g * A_tank * slosh_amp^2;

% Stokes boundary layer on side wall and floor
sh2       = sinh(k * h).^2;
int_wall  = 2 * pi * tank_r * (sinh(2 * k * h) / (4 * k) - h / 2) ./ sh2;
int_floor = A_tank ./ sh2;
P_wall    = 0.5 * prop_density * sqrt(nu * omega / 2) .* ...
    omega.^2 * slosh_amp^2 .* (int_wall + int_floor);
zeta_wall = P_wall ./ (2 * omega * E_modal);

% Ring baffles
d = h - z_baf;   % depth below quiescent surface (< 0 if above)

% vertical displacement amplitude across the ring
X = (sinh(k * min(z_baf, h)) ./ sinh(k * h)) .* eta_ring;

% fraction of the cycle each ring point is under the surface
s = 0.5 + asin(max(min(d ./ eta_ring, 1), -1)) / pi;

% X -> 0 at the floor and the J0 node (r = 0.63 R): C_D * X^3 -> 0
Cd  = 15 * (pi * max(X, eps) / baffle_w).^(-1/2);
eta = baffleThicknessFactor(baffle_t / baffle_w, 2 * X / baffle_w);

ring_int = trapz(r_ring(:), s .* eta .* Cd .* X.^3 .* 2 * pi .* r_ring, 3);
D_baffle = (4 * prop_density / 3) * omega.^2 .* sum(ring_int, 2);
zeta_baffle = D_baffle / (4 * pi * E_modal);

zeta_axial = reshape(zeta_wall + zeta_baffle, size(prop_mass_profile));
end

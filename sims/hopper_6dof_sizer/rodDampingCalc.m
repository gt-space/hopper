function [zeta_rods] = rodDampingCalc( ...
    prop_mass_profile, prop_density, tank_r, rods, slosh_amp)
% Lateral (m = 1) slosh damping from vertical members of the baffle mount
% (rods / flat bars running along the tank axis).
%
% Same energy method as the ring baffles (Miles 1958):
%   zeta = (energy dissipated per cycle) / (4*pi * modal energy)
%
% Mode shape (potential flow), k = 1.8412 / R:
%   eta(r, th) = eta_w * J1(k r) / J1(k R) * cos(th)
%   horizontal velocity amplitude at depth z (from the floor):
%     u_r  = omega * eta_w / (k J1(kR)) * k J1'(k r) cos(th) * cosh(k z) / sinh(k h)
%     u_th = omega * eta_w / (k J1(kR)) * J1(k r) / r  sin(th) * cosh(k z) / sinh(k h)
%   E = 1/2 rho g eta_w^2 * pi * Int_0^R (J1(kr)/J1(kR))^2 r dr
%
% Each member: Morison cross-flow drag per unit length 1/2 rho C_D D |u| u,
% dissipating (4/3) rho C_D D U^3 / omega per cycle per unit length, with
% the Keulegan-Carpenter drag coefficient used for the ring baffles,
%   C_D = 15 * (U T / D)^(-1/2)   (flat plates, KC small)
% capped at C_D_max for a bluff section. D is the width presented to the
% flow, which depends on the flow direction for a flat bar.
%
% rods: struct array, one entry per member, fields
%   r      radial position of the member centreline, m
%   th     angular position, rad (only relative angles matter)
%   width  section width, m      (round rod: diameter)
%   thick  section thickness, m  (round rod: diameter)
%   orient 'radial' (wide face along r), 'tangential', or 'round'
%   z0, z1 axial extent from the tank floor, m
% slosh_amp - wall slosh amplitude eta_w, m (default 0.1 * tank_r)

if nargin < 5, slosh_amp = 0.1 * tank_r; end

g       = 9.81;
k       = 1.8412 / tank_r;
A_tank  = pi * tank_r^2;
Cd_max  = 2.0;
J1R     = besselj(1, k * tank_r);

rr = linspace(0, tank_r, 400);
E_modal = 0.5 * prop_density * g * slosh_amp^2 * pi * ...
    trapz(rr, (besselj(1, k * rr) / J1R).^2 .* rr);

h_vec     = max((prop_mass_profile / prop_density) / A_tank, 1e-3);
zeta_rods = zeros(size(h_vec));

for i = 1:numel(h_vec)
    h     = h_vec(i);
    omega = sqrt(g * k * tanh(k * h));
    D_cyc = 0;
    for j = 1:numel(rods)
        rd = rods(j);
        z  = linspace(rd.z0, min(rd.z1, h), 60);
        if z(end) <= z(1), continue, end

        kr   = k * rd.r;
        dJ1  = 0.5 * (besselj(0, kr) - besselj(2, kr));
        J1_r = besselj(1, kr) / max(kr, eps);            % J1(kr) / (kr)
        amp  = slosh_amp * omega / J1R * cosh(k * z) / sinh(k * h);
        ur   = amp * abs(dJ1  * cos(rd.th));
        ut   = amp * abs(J1_r * sin(rd.th));

        % drag of each flow component on the width it sees
        switch rd.orient
            case 'radial',     Dr = rd.thick; Dt = rd.width;
            case 'tangential', Dr = rd.width; Dt = rd.thick;
            otherwise,         Dr = rd.width; Dt = rd.width;
        end
        D_cyc = D_cyc + segDissipation(ur, Dr, z, omega, prop_density, Cd_max) ...
                      + segDissipation(ut, Dt, z, omega, prop_density, Cd_max);
    end
    zeta_rods(i) = D_cyc / (4 * pi * E_modal);
end
end

function D = segDissipation(U, Dw, z, omega, rho, Cd_max)
% energy per cycle from quadratic drag along one member
UT = 2 * pi * max(U, eps) / omega;                % U_m * T
Cd = min(15 * (UT / Dw).^(-1/2), Cd_max);
D  = trapz(z, (4 / 3) * rho * Cd * Dw .* U.^3 / omega);
end

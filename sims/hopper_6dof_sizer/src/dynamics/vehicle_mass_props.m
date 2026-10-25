function [cg, I, dcg_dMox, dcg_dMfu, dI_dMox, dI_dMfu] = vehicle_mass_props(DRY, TG, m_ox, m_fu)
% Vehicle CG and inertia: dry block plus the two propellant columns.
%
%   DRY - dry block from mass_properties_setup:
%         mass [kg], cg [3x1 m, datum frame: x up the vehicle axis from
%         the base], I [3x3 kg m^2, about the dry CG, body axes x/y/z]
%   TG  - tank geometry: tank_r, ox_floor, fu_floor [m], ox_rho, fu_rho
%
%   cg  - vehicle CG [3x1 m], datum frame
%   I   - vehicle inertia about cg [3x3 kg m^2], body axes
%   dcg_dMox, dcg_dMfu - d(cg_x)/dm per tank [m/kg]
%   dI_dMox,  dI_dMfu  - dI/dm per tank [3x3 m^2]
%
% Each liquid sits as a flat column on its tank floor. The derivatives are
% numerical differences of this same function, so the integrated sim values
% and the directly evaluated ones cannot drift apart.

[cg, I] = props(DRY, TG, m_ox, m_fu);

if nargout > 2
    h = 1e-4;   % kg
    [dcg_dMox, dI_dMox] = numderiv(@(m) props(DRY, TG, m, m_fu), m_ox, h);
    [dcg_dMfu, dI_dMfu] = numderiv(@(m) props(DRY, TG, m_ox, m), m_fu, h);
end
end

function [cg, I] = props(DRY, TG, m_ox, m_fu)
R = TG.tank_r;
A = pi * R^2;
h_ox = m_ox / (TG.ox_rho * A);
h_fu = m_fu / (TG.fu_rho * A);
r_ox = [TG.ox_floor + h_ox / 2; 0; 0];
r_fu = [TG.fu_floor + h_fu / 2; 0; 0];

M  = DRY.mass + m_ox + m_fu;
cg = (DRY.mass * DRY.cg + m_ox * r_ox + m_fu * r_fu) / M;

I = DRY.I + parallel_axis(DRY.mass, DRY.cg - cg) ...
  + column(m_ox, R, h_ox) + parallel_axis(m_ox, r_ox - cg) ...
  + column(m_fu, R, h_fu) + parallel_axis(m_fu, r_fu - cg);
end

function J = column(m, R, h)
% solid cylinder about its own CG, axis along body x
J = diag([m * R^2 / 2, m * (3 * R^2 + h^2) / 12, m * (3 * R^2 + h^2) / 12]);
end

function J = parallel_axis(m, d)
J = m * ((d' * d) * eye(3) - d * d');
end

function [dcg, dI] = numderiv(f, m, h)
if m > h
    [cp, Ip] = f(m + h);
    [cm, Im] = f(m - h);
    dcg = (cp(1) - cm(1)) / (2 * h);
    dI  = (Ip - Im) / (2 * h);
else
    % one-sided at an empty tank
    [c0, I0] = f(m);
    [c1, I1] = f(m + h);
    [c2, I2] = f(m + 2 * h);
    dcg = (-3 * c0(1) + 4 * c1(1) - c2(1)) / (2 * h);
    dI  = (-3 * I0 + 4 * I1 - I2) / (2 * h);
end
end

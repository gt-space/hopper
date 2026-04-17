function [zeta_axial] = axialDampingCalc( ...
    prop_mass_profile, prop_density, nu, tank_r, baffle_w, Nb, tank_h)

g      = 9.81;
A_tank = pi * tank_r^2;
r_gap  = tank_r - baffle_w;
A_gap  = pi * r_gap^2;
sigma  = 1 - (A_gap / A_tank);
Cd     = 1.5;

h_vec      = (prop_mass_profile / prop_density) / A_tank;
zeta_axial = zeros(size(h_vec));
z_baf      = linspace(0, tank_h, Nb);

baffle_spacing = tank_h / max(Nb - 1, 1);

for i = 1:length(h_vec)
    h = h_vec(i);

    delta = 0.2 * baffle_spacing;
    
    omega = sqrt(g / max(h, 1e-6));
   
    alpha_part = 0.5 * (1 - exp(-2 * h / tank_h));   % ramps up with fill,
                                                       
    m_modal = prop_density * A_tank * h * alpha_part;
  
    E_modal = 0.5 * m_modal * omega^2;

    % Smooth-wall axial damping (viscous BL on cylindrical wall) 
    % SP-106 axial viscous term
    zeta_smooth = (1 / (2 * tank_r)) * sqrt(nu / (2 * omega));

    % Axial mode shape 
    % First axial mode: sinusoidal pressure variation, 
    % velocity amplitude proportional to sin(pi*z/h)
    % Maximum at surface (z=h), zero at bottom
    phi = sin((pi / 2) * (z_baf / max(h, 1e-6)));
    phi = min(phi, 1);   % cap at free surface value

   
    w = 1 ./ (1 + exp((z_baf - h) / delta));

    zeta_baffle = 0;
    for j = 1:Nb
        phi_sq = phi(j)^2;
        weight = phi_sq / (1 + phi_sq);

        dE_i       = w(j) * Cd * sigma^1.2 * weight * E_modal;
        zeta_baffle = zeta_baffle + dE_i / (4 * pi * E_modal);
    end

    zeta_axial(i) = zeta_smooth + zeta_baffle;
end
end
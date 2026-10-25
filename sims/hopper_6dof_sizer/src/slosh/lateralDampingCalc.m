function [damping_ratios] = lateralDampingCalc( ...
    prop_mass_profile, prop_density, nu, tank_r, baffle_w, Nb, tank_h, ...
    baffle_t, slosh_amp)

% baffle_t  - baffle thickness [m] (default 0 -> thin plate, no correction)
% slosh_amp - free surface slosh amplitude at tank wall, y_s [m]
%             (default 0.1 * tank_r)
if nargin < 8, baffle_t  = 0;            end
if nargin < 9, slosh_amp = 0.1 * tank_r; end

g = 9.81;
eig_1 = 1.841; % 1st mode of oscillation for cylindrical tanks

prop_heights = (prop_mass_profile / prop_density) / (pi * tank_r^2);

% smooth-wall correlation
alpha_lut = [0 0.25 0.5 1 1.5 2 10];
c_1_lut   = [5 4.5 2.5 1.8 1.75 1.75 1.75];

damping_ratios = zeros(size(prop_heights));

z_baf = linspace(0, tank_h, Nb);

% geometry constants
A_tank = pi * tank_r^2;
r_gap = tank_r - baffle_w;
A_gap = pi * r_gap^2;
sigma = 1 - (A_gap / A_tank);

k = eig_1 / tank_r;

Cd = 1.7;

for i = 1:length(prop_heights)

    h = prop_heights(i);

    % Slosh frequency (SP-106)
   
    omega = sqrt((g * eig_1 / tank_r) * tanh(eig_1 * h / tank_r));

    % Smooth-wall damping
  
    alpha_fill = h / tank_r;

    c_1 = interp1(alpha_lut, c_1_lut, alpha_fill, 'linear', 'extrap');

    zeta_smooth = (c_1 / (2 * tank_r)) * sqrt(nu / (2 * omega));

    % potential-flow mode shape at each ring, 1 at the free surface,
    % 1/cosh(k h) at the floor (rings above the liquid see the surface value)
    phi = cosh(k * min(z_baf, h)) ./ cosh(k * h);

    alpha_part = 0.5 * (1 - exp(-2*h/tank_r));

    m_modal = prop_density * A_tank * h * alpha_part;

    % modal energy

    E_modal = 0.5 * m_modal * omega^2;
  
    delta = max(0.03 * tank_h, 0.05 * h);

    w = 1 ./ (1 + exp((z_baf - h) / delta));

    % Baffle thickness effectiveness (NASA TN D-3716)
    % mean double amplitude at baffle edge around the ring:
    % A = (4/pi) * y_s * exp(-1.84 d/a), d = baffle depth below surface
    d_baf = max(h - z_baf, 0);
    A_edge = (4 / pi) * slosh_amp * exp(-eig_1 * d_baf / tank_r);
    eta = baffleThicknessFactor(baffle_t / baffle_w, A_edge / baffle_w);

    zeta_energy = 0;

    for j = 1:Nb

        % mode participation (NO normalization)
        phi_sq = phi(j)^2;

        % local participation weight (relative, not normalized)
        weight = phi_sq / (1 + phi_sq);  

        % dissipation as FRACTION of modal energy
        dE_i = eta(j) * w(j) * Cd * sigma^1.2 * weight * E_modal;

        zeta_energy = zeta_energy + dE_i;

    end

  
    zeta_baffle = zeta_energy / (4 * pi * E_modal);

    % total damping
    
    damping_ratios(i) = zeta_smooth + zeta_baffle;

    % disp(zeta_baffle)

end

end
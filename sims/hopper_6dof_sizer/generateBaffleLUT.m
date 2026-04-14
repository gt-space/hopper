
ox_density = 1141;   % kg/m^3
fu_density = 800;    % kg/m^3

fu_nu = 1.3e-6;  % m^2/s FU
ox_nu = 0.17e-6; % m^2/s OX

tank_r = TANKS.singular.radius;
ox_baffle_width = TANKS.singular.radius * 0.35;
fu_baffle_width = TANKS.singular.radius * 0.35;
ox_baffle_number = 10;
fu_baffle_number = 10;

ox_tank_h = TANKS.singular.oxidizer.h;
fu_tank_h = TANKS.singular.fuel.h;

ox_mass_profile = linspace(1, IN.propulsion.oxidizer_mass, 200);
fu_mass_profile = linspace(1, IN.propulsion.fuel_mass, 200);

[ox_damping_ratios] = lateralDampingCalc(ox_mass_profile, ox_density, ox_nu, tank_r, ox_baffle_width, ox_baffle_number, ox_tank_h);

[fu_damping_ratios] = lateralDampingCalc(fu_mass_profile, fu_density, fu_nu, tank_r, fu_baffle_width, fu_baffle_number, fu_tank_h);

figure
plot(ox_mass_profile, ox_damping_ratios)
xlabel("Ox mass")
ylabel("Lateral Damping Ratio")
title("Lateral Damping Ratio over Ox mass, Number of Baffles: " + ox_baffle_number + " , Baffle width (m): " + ox_baffle_width)

figure
plot(fu_mass_profile, fu_damping_ratios)
xlabel("Fuel mass")
ylabel("Lateral Damping Ratio")
title("Lateral Damping Ratio over Fuel mass, Number of Baffles: " + fu_baffle_number + " , Baffle width(m): " + fu_baffle_width)
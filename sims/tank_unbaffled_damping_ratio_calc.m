
liq_depths = linspace(0.01,0.5463,200); % m
g = 9.81; % m/s^2
eig_1 = 1.841;
nu = 1.3e-6;% m^2/s FU
% nu = 0.17e-6; % m^2/s OX
tank_r = 0.1206;
tank_h = 0.4687;

alpha_lut = [0 0.25 0.5 1 1.5 2 10];
c_1_lut = [5 4.5 2.5 1.8 1.75 1.75 1.75];

damping_ratio = size(liq_depths);

for i = 1:length(liq_depths)
    omega = sqrt(g * eig_1 * tanh(eig_1*liq_depths(i)));
    alpha = liq_depths(i) / tank_r;
    c_1 = interp1(alpha_lut, c_1_lut, alpha);
    damping_ratio(i) = 1/(2*tank_r) * sqrt(nu/(2*omega)) * c_1;
end

plot(liq_depths, damping_ratio)
xlabel("Propellant Height (m)")
ylabel("Damping Ratio")
title("Damping Ratio over Propellant Height")
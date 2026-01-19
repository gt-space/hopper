% Function Testing
thrust_target = 200 * 4.44822; % lbf --> N
ox_name = "N2O";
fu_name = "Ethanol";
fu_temp = 293; % K
ox_vap_P = 500 * 6894.76; % psia --> Pa
ox_tank_P = 500 * 6894.76; % psia --> Pa
fu_tank_P = 500 * 6894.76; % psia --> Pa
ox_sys_CdA = 2.0975E-5; % m^2
fu_sys_CdA = 7.588E-6; % m^2
eta_cstar = 0.85;
At = 1.150 * 0.00064516; % in^2 --> m^2
eps = 3.356;

[ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR, ox_valve_CdA, fu_valve_CdA] = prop_system(thrust_target, ox_name, fu_name, fu_temp, ox_vap_P, ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps);

Pc = Pc / 6895.76; % Pa --> psia

% Prints
fprintf('OX Mdot : %.4f kg/s\n', ox_mdot)
fprintf('FU Mdot : %.4f kg/s\n', fu_mdot)
fprintf('Total Mdot : %.4f kg/s\n', tot_mdot)
fprintf('Pc : %.3f psia\n', Pc)
fprintf('Thrust : %.2f lbf\n', thrust / 4.44822)
fprintf('MR : %.2f\n', MR)
fprintf('OX Valve CdA : %.10f m^2\n', ox_valve_CdA)
fprintf('FU Valve CdA : %.10f m^2\n', fu_valve_CdA)

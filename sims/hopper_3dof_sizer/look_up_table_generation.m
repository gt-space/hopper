thrust_target = 200:540 * 4.44822; % lbf --> N
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

mdot_data = cell(length(thrust_target),2);

for i = 1:length(thrust_target)
    [ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR, ox_valve_CdA, fu_valve_CdA] = prop_system(thrust_target(i), ox_name, fu_name, fu_temp, ox_vap_P, ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps);
    mdot_data{i,1} = thrust_target(i);
    mdot_data{i,2} = tot_mdot;
end

writecell(mdot_data, 'mdot_lookup.xlsx')

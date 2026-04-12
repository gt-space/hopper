thrust_target = 200:540 * 4.44822; % lbf --> N
ox_name = "O2";
fu_name = "Dodecane";
fu_temp = 293; % K
ox_temp = 90; % K
ox_tank_P = 500 * 6894.76; % psia --> Pa
fu_tank_P = 500 * 6894.76; % psia --> Pa
ox_sys_CdA = 1.7827-5; % m^2
fu_sys_CdA = 1.1E-5; % m^2
eta_cstar = 0.85;
At = 1.056 * 0.00064516; % in^2 --> m^2
eps = 2.93;

mdot_data = cell(length(thrust_target),3);

for i = 1:length(thrust_target)
    [ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR, ox_valve_CdA, fu_valve_CdA] = prop_system(thrust_target(i), ox_name, fu_name, fu_temp, ox_temp, ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps);
    mdot_data{i,1} = thrust_target(i);
    mdot_data{i,2} = ox_mdot;
    mdot_data{i,3} = fu_mdot;
end

writecell(mdot_data, 'mdot_lookup.xlsx');

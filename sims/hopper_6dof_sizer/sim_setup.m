thrust_target = 200 * 4.44822; % lbf --> N
ox_name = "O2";
fu_name = "Dodecane";
fu_temp = 293; % K
ox_temp = 90; % K
ox_tank_P = 500 * 6894.76; % psia --> Pa
fu_tank_P = 500 * 6894.76; % psia --> Pa
ox_sys_CdA = 1.7827-5; % m^2
fu_sys_CdA = 1.1E-5; % m^2
eta_cstar = 0.85;
At = 1.138 * 0.00064516; % in^2 --> m^2
eps = 3.655;

main();

load_lookup();

load("cg_I_LUT.mat");

cg_moi_test();
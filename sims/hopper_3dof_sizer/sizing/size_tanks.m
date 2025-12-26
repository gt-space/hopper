function tanks = size_tanks(IN, prop)

%% Densities
rho_ox = 750;  % kg/m^3 (subcooled N2O)
rho_fu = 785;  % kg/m^3 (IPA)

V_ox = prop.m_ox / rho_ox;
V_fu = prop.m_fu / rho_fu;

%% Pressure assumptions
Pc = IN.propulsion.chamber_pressure;
P_tank = Pc * (1 + IN.propulsion.injector_dp_frac);

%% Material
sigma = 2.5e8; % Pa
rho_mat = 2800;
FS = IN.const.FS_struct;

%% Cylindrical tank proxy
function m = tank_mass(V)
    r = (V/pi)^(1/3);
    L = V / (pi*r^2);
    t = P_tank * r / (sigma / FS);
    m = 2*pi*r*L*t*rho_mat;
end

m_ox_tank = tank_mass(V_ox);
m_fu_tank = tank_mass(V_fu);

%% COPV sizing
if strcmp(IN.press.mode,'copv')
    P0 = IN.press.copv.max_pressure;
    Vcopv = (P_tank * (V_ox+V_fu)) / P0;
else
    Vcopv = 0;
end

%% Pack output
tanks = struct();
tanks.ox.volume = V_ox;
tanks.fu.volume = V_fu;
tanks.ox.mass = m_ox_tank;
tanks.fu.mass = m_fu_tank;
tanks.copv.volume = Vcopv;

end

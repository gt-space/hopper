function TANKS = size_tanks(IN)
% Sizes propellant tanks for multiple geometries
% Outputs all options: singular, clustered, concentric

%% Material properties (Al6061-T6)
mat.rho = 2700;              % kg/m^3
mat.sigma_y = 276e6;         % Pa

FS = IN.const.FS_struct;
P = 5.9e6;                   % Pa, assumed MEOP (850 psi)

%% Propellant properties
rho_ox = 750;    % kg/m^3 (N2O  intial) % TODO update from fill sim
rho_fu = 784;    % kg/m^3 (IPA)

m_ox = IN.propulsion.oxidizer_mass;
m_fu = IN.propulsion.fuel_mass;

V_ox = m_ox / rho_ox;
V_fu = m_fu / rho_fu;

TANKS.volumes.oxidizer = V_ox;
TANKS.volumes.fuel     = V_fu;

%% Helper function
tank_mass = @(V, r) ...
    (2*pi*r^2 + 2*pi*r*(V/(pi*r^2))) * ...
    (P*r/(mat.sigma_y/FS)) * mat.rho;

% TODO table of all possible pipe sizes that we can uses
%% ---------- 1. Singular (stacked) ----------
r = IN.tanks.max_radius;

m_ox_sing = tank_mass(V_ox, r);
m_fu_sing = tank_mass(V_fu, r);

TANKS.singular.oxidizer.mass = m_ox_sing;
TANKS.singular.fuel.mass     = m_fu_sing;
TANKS.singular.total_mass    = m_ox_sing + m_fu_sing;

%% ---------- 2. Clustered ----------
N = 3; % TODO maybe make this an input or size it
r_c = r / sqrt(N); % fit tanks within max radius - maybe chance to account for totally different hopper shape

m_ox_clust = N * tank_mass(V_ox/N, r_c);
m_fu_clust = N * tank_mass(V_fu/N, r_c);

TANKS.clustered.oxidizer.mass = m_ox_clust;
TANKS.clustered.fuel.mass     = m_fu_clust;
TANKS.clustered.total_mass    = m_ox_clust + m_fu_clust;

%% ---------- 3. Concentric ----------
r_inner = 0.6 * r;
r_outer = r;

m_fu_con = tank_mass(V_fu, r_inner);
m_ox_con = tank_mass(V_ox, r_outer);

TANKS.concentric.oxidizer.mass = m_ox_con;
TANKS.concentric.fuel.mass     = m_fu_con;
TANKS.concentric.total_mass    = m_ox_con + m_fu_con;

end

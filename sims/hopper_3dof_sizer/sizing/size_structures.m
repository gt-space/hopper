function structures = size_structures(IN, vehicle)

m = vehicle.mass_wet;
h = IN.mission.target_altitude;

%% Landing legs (Euler buckling proxy)
E = 69e9; % Pa (Al)
L = 0.5;  % m
I = 1e-8; % m^4

Pcr = pi^2 * E * I / L^2;
assert(Pcr > m*IN.const.g0, 'Landing leg buckling');

m_legs = 0.05 * m * IN.legs.tip_factor;

%% TVC actuators
m_tvc = 0.02 * m;

%% Pack
structures = struct();
structures.legs = m_legs;
structures.tvc = m_tvc;
structures.total = m_legs + m_tvc;

end

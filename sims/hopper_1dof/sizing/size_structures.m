function [STRUCT, VEH] = size_structures(IN, VEH)
% Sizes landing legs and intertank structure

%% Material Al6061
rho = 2700;        % kg/m^3
E   = 69e9;        % Pa

g = IN.const.g0;
FS = IN.const.FS_struct;

%% Vehicle loads
W = VEH.mass.wet * g * 2;
H = 0.81;   % engine bay height [m] (adjust if needed)
leg_length = H / cosd(30);

%% Leg sizing
N_legs = 4;
F_leg = W / N_legs * IN.legs.tip_factor;

% Buckling-limited rod
L = leg_length;
r = 0.001; % m (guess)
I = pi*r^4/4;

Pcr = pi^2 * E * I / L^2;

while Pcr < FS * F_leg
    r = r * 1.05;
    I = pi*r^4/4;
    Pcr = pi^2 * E * I / L^2;
end

A = pi*r^2;
m_leg = A * L * rho;

STRUCT.legs.mass = N_legs * m_leg * 1.5;

%% Intertank structure
STRUCT.intertank.mass = 3; % guess

STRUCT.total = STRUCT.legs.mass + STRUCT.intertank.mass;

VEH.mass.dry = VEH.mass.dry + STRUCT.total;

end

% IMPORTANT: This code used CoolProp and as such you must have a CoolProp
% and a compatible version of Python downloaded and added to path for
% MATLAB to properly. 

function [CoM,MoI,mass_IPA,mass_N2O,H_N2O] = HopperCoMMoI(Thrust,mass_IPA,mass_N2O,H_N2O)
% assuming thrust required given in lbf
Thrust = Thrust; %lbf
% converting to Newtons for simplicity in calculations
T = Thrust*4.4482216153; %N
% masses are in kg
% Volumes are in m^3
% enthalpy are in J/kg
% initialize a note cell array with notes for future actions needed in code
note = [];

%% Thrust to mdots (the bad way)
% % specific impulse of our IPA/N2O engine in seconds
% Isp = 189; %sec
% note = [note;{'confirm Isp'}];
% % sea level gravity
% g0 = 9.81; %m/s^2
% % total mass flow rate from ISP and thrust
% mdot = T/(Isp*g0); %kg/s
% % mixture ratio of propellants by mass
% MR = 3.8; %by mass  %3.8 is for 0.57:0.15
% note = [note;{'confirm MR'}];
% % mass flow rated of each propellant from total mass flow rate and mixture ratio
% mdot_N2O = MR/(MR+1)*mdot; %kg/s
% mdot_IPA = 1/(MR+1)*mdot; %kg/s

%% Thrust to mdots (the good way)
% tank to injector and injector to chamber CdA's
CdA_fs = []
CdA_inj = []
% fluid properties
rho = []
gamma = []
MW = []
R = 8.314/MW; %J/kg/K
% Tank ullage pressure (N2O gas pressure) from total N2O specific enthalpy and quality
    rho_IPA = 786; %kg/m^3
    note = [note;{'confirm IPA density'}];
        % find volume of IPA from mass and density
    V_IPA = mass_IPA/rho_IPA; %m^3
        % total volume of both tanks
    V_tot = 62e-3; %m^3
        % find volume of N2O from total volume minus IPA volume
    V_N2O = V_tot-V_IPA; %m^3
        % find the overall density of the N2O liquid/gas mixture
    rho_N2O_tot = mass_N2O/V_N2O; %kg/m^3
    h_N2O_tot = H_N2O/mass_N2O; %J/kg
    quality = py.CoolProp.CoolProp.PropsSI('Q','H',h_N2O_tot,'D',rho_N2O_tot,'N2O');
p_tank = py.CoolProp.CoolProp.PropsSI('P','Q',quality,'H',h_N2O_tot,'N2O');
% rocket engine properties
A_t = []
A_e = []
Astar = A_t;
cstar = 0.85;
syms mdot p_inj p_c M_e T_e p_e u_e T_c
% CdA eqn for feed system
eqn1 = CdA_fs == mdot/sqrt(2*rho*(p_tank-p_inj));
% CdA eqn for injector
eqn2 = CdA_inj == mdot/sqrt(2*rho*(p_inj-p_c));
% cstar eqn
eqn3 = mdot == p_c*Astar/cstar;
% Rocket thrust eqn
eqn4 = T == mdot*u_e + (p_e-p_c)*A_e;
% exit velocity eqn
eqn5 = u_e == M_e*sqrt(gamma*R*T_e);
% Isentropic area ratio eqn
eqn6 = A_e/Astar == ((gamma+1)/2)^(-(gamma+1)/(2*gamma-2))*((1+(gamma-1)/2*M_e^2)^((gamma+1)/(2*gamma-2))/M_e);
% isentropic pressure ratio eqn
eqn7 = p_c == p_e*(1+(gamma-1)/2*M_e^2)^(gamma/(gamma-1));

eqns = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7];
values = solve(eqns, [mdot p_inj p_c M_e T_e p_e u_e]);
clear mdot
mdot = values.mdot;

MR = 3.8; %by mass  %3.8 is for 0.57:0.15
note = [note;{'confirm MR'}];
% mass flow rated of each propellant from total mass flow rate and mixture ratio
mdot_N2O = MR/(MR+1)*mdot; %kg/s
mdot_IPA = 1/(MR+1)*mdot; %kg/s

%% Mdots to Updated Mass
% time step used in simulation
dt = 20e-5; %sec
note = [note;{'confirm dt'}];
% updating the mass of each propellant
mass_N2O = mass_N2O - mdot_N2O*dt; %kg
mass_IPA = mass_IPA - mdot_IPA*dt; %kg

%% Mass to Thermodynamic State

% intial volumes
Vi_IPA = 12e-3; %m^3
Vi_N2O = 50e-3; %m^3
note = [note;{'can the initial volumes be set outside of the code?'}];
note = [note;{'confirm Vi/V_tot'}];

% total volume of both tanks
V_tot = 62e-3; %m^3

rho_IPA = 786; %kg/m^3
note = [note;{'confirm IPA density'}];

% find volume of IPA from mass and density
V_IPA = mass_IPA/rho_IPA; %m^3
% find volume of N2O from total volume minus IPA volume
V_N2O = V_tot-V_IPA; %m^3
% find the overall density of the N2O liquid/gas mixture
rho_N2O_tot = mass_N2O/V_N2O; %kg/m^3

% initial masses
mi_IPA = 9; %kg
mi_N2O = 40; %kg
note = [note;{'can the initial masses be set outside of the code?'}];

% enthalpy of the N2O liquid from mixture density and liquid quality
h_N2O_liquid = py.CoolProp.CoolProp.PropsSI('H','D',rho_N2O_tot,'Q',0,'N2O'); %J/kg
note = [note;{'I think this step may be wrong'}];
% updating N2O mixture enthalpy
H_N2O = H_N2O-h_N2O_liquid*mdot_N2O*dt;
% updating N2O mixture specific enthalpy
h_N2O_tot = H_N2O/mass_N2O;

% quality of N2O mixture
quality = py.CoolProp.CoolProp.PropsSI('Q','H',h_N2O_tot,'D',rho_N2O_tot,'N2O');
% finding mass of N2O liquid and gas from quality and total mass
mass_N2O_gas = quality*mass_N2O;
mass_N2O_liq = (1-quality)*mass_N2O;

% N2O gas pressure from quality and mixture specific enthalpy
P_N2O = py.CoolProp.CoolProp.PropsSI('P','Q',quality,'H',h_N2O_tot,'N2O');
% density of N2O liquid and gas from N2O pressure and l/g qualities
rho_N2O_liq = py.CoolProp.CoolProp.PropsSI('D','Q',0,'P',P_N2O,'N2O');
rho_N2O_gas = py.CoolProp.CoolProp.PropsSI('D','Q',1,'P',P_N2O,'N2O');
note = [note;{'I think this step may be wrong'}];
% volumes of liquid and gas N2O from mass and densities
V_N2O_liq = mass_N2O_liq/rho_N2O_liq; %m^3
V_N2O_gas = mass_N2O_gas/rho_N2O_gas; %m^3

%% Thermodynamic State to Center of Mass
% center of mass relative to bottom of tanks
% assuming the vehicle is symetric about the z-axis

% tank dimensions ignoring thickness and end caps
r = 2; %inch
R = 4; %inch
H_tank = 1.9; %m
h_tank = 0.85*H_tank; %m
note = [note;{'update with final dimensions and include endcaps'}];
% convert everything to meters
r = 0.0254*r; %m
R = 0.0254*R; %m

% area of tanks
a = pi*r^2; %m^2
A = pi*R^2 - a; %m^2
top_A = pi*R^2; %m^2

% height of liquids (assumes liquid N2O never reach top of inner tank
height_IPA = V_IPA/a; %m
height_N2O_liq = V_N2O_liq/A; %m
height_N2O_gas_ring = h_tank-height_N2O_liq; %m
height_N2O_gas_top = H_tank-h_tank; %m
height_N2O_gas_inner = h_tank-height_IPA; %m

% center of mass of individual components relative to bottom of tanks
% ipa
CoM_IPA = [0;0;height_IPA/2]; %m
% N2O liquid
CoM_N2O_liq = [0;0;height_N2O_liq/2]; %m
% N2O gas in outer tank below height of inner tank
CoM_N2O_gas_ring = [0;0;height_N2O_liq + height_N2O_gas_ring/2]; %m
% N2O gas in outer tank above height of inner tank
CoM_N2O_gas_top = [0;0;h_tank + height_N2O_gas_top/2]; %m
% N2O gas in inner tank
CoM_N2O_gas_inner = [0;0;height_IPA + height_N2O_gas_inner/2]; %m
% inner tank plunger
CoM_plunger = [0;0;height_IPA]; %m
% the center of mass of the entire structure
CoM_structure = [0;0;-0.15]; %m
note = [note;{'update when the final CAD is done. For now it just ~6 inches below tanks'}];

% masses of individual compenents in kilograms
%mass_IPA; %kg
%mass_N2O_liq; %kg
% volumes of N2O gas
V_N2O_gas_ring = A*(h_tank-height_N2O_liq); %m^3
V_N2O_gas_top = top_A*(H-h); %m^3
V_N2O_gas_inner = a*(h_tank-height_IPA); %m^3
% masses of N2O gas
mass_N2O_gas_ring = rho_N2O_gas*V_N2O_gas_ring; %kg
mass_N2O_gas_top = rho_N2O_gas*V_N2O_gas_top; %kg
mass_N2O_gas_inner = rho_N2O_gas*V_N2O_gas_inner; %kg
% mass of plunger
mass_plunger = 0.0225; %kg
note = [note;{'update with real mass'}];
% mass of total structure excluding plunger
mass_structure = 15; %kg
note = [note;{'update with real mass'}];

% Total Center of Mass
CoM_top = mass_IPA*CoM_IPA + mass_N2O_liquid*CoM_N2O_liq + mass_N2O_gas_ring*CoM_N2O_gas_ring + mass_N2O_gas_top*CoM_N2O_gas_top + mass_N2O_gas_inner*CoM_N2O_gas_inner + mass_plunger*CoM_plunger + mass_structure*CoM_structure; %kg*m
CoM_bot = mass_IPA + mass_N2O_liquid + mass_N2O_gas_ring + mass_N2O_gas_top + mass_N2O_gas_inner + mass_plunger + mass_structure; %kg
CoM = CoM_top/CoM_bot; %m

% distances from center of mass of each component center of mass
d_IPA = CoM_IPA - CoM; %m
d_N2O_liq = CoM_N2O_liq - CoM; %m
d_N2O_gas_ring = CoM_N2O_gas_ring - CoM; %m
d_N2O_gas_top = CoM_N2O_gas_top - CoM; %m
d_N2O_gas_inner = CoM_N2O_gas_inner - CoM; %m
d_plunger = CoM_plunger - CoM; %m
d_structure = CoM_structure - CoM; %m

%% Thermodynamic State to Moment of Inertia

% moments of inertia of each component about their centers of mass
% assuming solid cylinders and thick walled cylindrical tubes for fluids and thin disk for plunger
% x and y are in the horizontal directions and equivalent due to symmetry, z in the up and down direction
% MoI matrices from https://en.wikipedia.org/wiki/List_of_moments_of_inertia
MoI_IPA = diag([1/12*mass_IPA*(3*r^2+height_IPA^2),1/12*mass_IPA*(3*r^2+height_IPA^2),1/2*mass_IPA*r^2]); %kg*m^2
MoI_N2O_liq = diag([1/12*mass_N2O_liq*(3*(r^2+R^2)+height_N2O_liq^2),1/12*mass_N2O_liq*(3*(r^2+R^2)+height_N2O_liq^2),1/2*mass_N2O_liq*(r^2+R^2)]); %kg*m^2
MoI_N2O_gas_ring = diag([1/12*mass_N2O_gas_ring*(3*(r^2+R^2)+height_N2O_gas_ring^2),1/12*mass_N2O_gas_ring*(3*(r^2+R^2)+height_N2O_gas_ring^2),1/2*mass_N2O_gas_ring*(r^2+R^2)]); %kg*m^2
MoI_N2O_gas_top = diag([1/12*mass_N2O_gas_top*(3*R^2+height_N2O_gas_top^2),1/12*mass_N2O_gas_top*(3*R^2+height_N2O_gas_top^2),1/2*mass_N2O_gas_top*R^2]); %kg*m^2
MoI_N2O_gas_inner = diag([1/12*mass_N2O_gas_inner*(3*r^2+height_N2O_gas_inner^2),1/12*mass_N2O_gas_inner*(3*r^2+height_N2O_gas_inner^2),1/2*mass_N2O_gas_inner*r^2]); %kg*m^2
MoI_plunger = diag([1/4*mass_plunger*r^2,1/4*mass_plunger*r^2,1/2*mass_plunger*r^2]); %kg*m^2
MoI_structure = diag(['idk','have absolutely','no clue']); %kg*m^2
note = [note;{'update with moment of inertia when get CAD'}];

% moments of inertia of each component about the overall center of mass using paralled axis theorem
MoI_IPA_CoM = MoI_IPA + mass_IPA.*d_IPA.^2; %kg*m^2
MoI_N2O_liq_CoM = MoI_N2O_liq + mass_N2O_liq.*d_N2O_liq.^2; %kg*m^2
MoI_N2O_gas_ring_CoM = MoI_N2O_gas_ring + mass_N2O_gas_ring.*d_N2O_gas_ring.^2; %kg*m^2
MoI_N2O_gas_top_CoM = MoI_N2O_gas_top + mass_N2O_gas_top.*d_N2O_gas_top.^2; %kg*m^2
MoI_N2O_gas_inner_CoM = MoI_N2O_gas_inner + mass_N2O_gas_inner.*d_N2O_gas_inner.^2; %kg*m^2
MoI_plunger_CoM = MoI_plunger + mass_plunger.*d_plunger.^2; %kg*m^2
MoI_structure_CoM = MoI_structure + mass_structure.*d_structure.^2; %kg*m^2

% overall MoI from adding MoI of each component
MoI = MoI_IPA_CoM + MoI_N2O_liq_CoM + MoI_N2O_gas_ring_CoM + MoI_N2O_gas_top_CoM + MoI_N2O_gas_inner_CoM + MoI_plunger_CoM + MoI_structure_CoM; %kg*m^2

end

Cstar_actual = 
Cstar_ideal = 1502; %m/s
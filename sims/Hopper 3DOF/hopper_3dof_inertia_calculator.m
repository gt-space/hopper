function [Iyy, Iyy_dot, cg, mdot_n2o, mdot_ipa, m] = hopper_3dof_inertia_calculator(ipa_mass, n2o_mass,mdot, vehicle)

% Calculate Total Mass
m = vehicle.mass_d + ipa_mass + n2o_mass;

% Calculate Height of Propellants
n2o_height = vehicle.n2o_tank_height * (n2o_mass / vehicle.n2o_mass_0);
ipa_height = vehicle.ipa_tank_height * (ipa_mass / vehicle.ipa_mass_0);

% Calculate Propellant CGs
n2o_cg = vehicle.n2o_tank_bottom + n2o_height / 2;
ipa_cg = vehicle.ipa_tank_bottom + ipa_height / 2;

% Calculate Total CG
cg = (vehicle.mass_d * vehicle.cg_dry ...
    + n2o_mass * n2o_cg ...
    + ipa_mass * ipa_cg) / m;

% Calculate Propellant Mass Flow Rates
mdot_n2o = vehicle.mr * mdot / (1+vehicle.mr);
mdot_ipa = mdot / (1+vehicle.mr);

% Calculate Rate of Change of Center of Gravity
cg_dot = ( ...
    mdot_n2o*(n2o_cg - vehicle.cg_dry) + ...
    mdot_ipa*(ipa_cg -  vehicle.cg_dry)) / m;


% Calculate Moment of Inertia
Iyy = vehicle.Iyy_dry ...
    + vehicle.mass_d*(vehicle.cg_dry - cg)^2 ...
    + n2o_mass*(n2o_cg - cg)^2 ...
    + ipa_mass*(ipa_cg - cg)^2;

% Calculate Rate of Change of Moment of Inertia
Iyy_dot = mdot_n2o*(n2o_cg - cg)^2 ...
    + mdot_ipa*(ipa_cg - cg)^2 ...
    - 2 * (n2o_mass*(n2o_cg - cg) + ipa_mass*(ipa_cg - cg)) * cg_dot;

end
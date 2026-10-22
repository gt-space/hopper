% Vehicle mass properties for the sim: dry block (main.m dry-mass override)
% plus the propellant columns. See mass_properties_setup and
% MASS_PROPERTIES.md.

% Variables used directly by Simulink blocks and other scripts
ox_mass = IN.propulsion.oxidizer_mass;
fu_mass = IN.propulsion.fuel_mass;
engine_mass = IN.engine_mass;
engine_cg = IN.engine_cg;
mount_cg = IN.mount_cg;
tank_r = TANKS.singular.radius;

MP = mass_properties_setup(IN, TANKS, OUT.Vehicle.DryMass);

cg_init  = MP.cg_init;
MoI_init = MP.MoI_init;

% Derivative tables for CG // MOI Script v3 (replace the cg_I_LUT.mat
% values loaded by load_lookup, so they always match the current masses)
ox_masses      = MP.ox_masses;
fu_masses      = MP.fu_masses;
dcg_dMox_table = MP.dcg_dMox_table;
dcg_dMfu_table = MP.dcg_dMfu_table;
dI_dMox_table  = MP.dI_dMox_table;
dI_dMfu_table  = MP.dI_dMfu_table;

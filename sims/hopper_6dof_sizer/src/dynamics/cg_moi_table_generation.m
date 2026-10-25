% Write the CG / inertia derivative tables to cg_I_LUT.mat.
% Run after main (needs IN, TANKS, OUT). sim_setup and mc_sim_setup rebuild
% these tables themselves; this file keeps cg_I_LUT.mat in step for
% anything that loads it directly.

MP = mass_properties_setup(IN, TANKS, OUT.Vehicle.DryMass);

ox_masses      = MP.ox_masses;
fu_masses      = MP.fu_masses;
dcg_dMox_table = MP.dcg_dMox_table;
dcg_dMfu_table = MP.dcg_dMfu_table;
dI_dMox_table  = MP.dI_dMox_table;
dI_dMfu_table  = MP.dI_dMfu_table;

save(hopperRoot('data', 'cg_I_LUT.mat'), ...
    "ox_masses", "fu_masses", ...
    "dcg_dMox_table", "dcg_dMfu_table", ...
    "dI_dMox_table", "dI_dMfu_table");

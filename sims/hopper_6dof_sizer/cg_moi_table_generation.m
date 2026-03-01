
ox_mass = IN.propulsion.oxidizer_mass;
fu_mass = IN.propulsion.fuel_mass;

engine_mass = ENGINE.mass;
engine_cg = 0.3556;
ox_tank_mass = TANKS.singular.oxidizer.mass;
ox_tank_cg = 1.854;
ox_tank_wall_thick = 0.003175;
ox_tank_h = TANKS.singular.oxidizer.h;
fu_tank_mass = TANKS.singular.fuel.mass;
fu_tank_cg = 1.0668;
fu_tank_wall_thick = 0.003175;
fu_tank_h = TANKS.singular.fuel.h;
structures_mass =  OUT.Vehicle.MassDistribution.Structures;
structures_cg = 1.193;
avi_mass = AVI.mass;
avi_cg = 1.4224;
fluids_mass = OUT.VehicleFluids.FinalMass;
fluids_cg = 1.0668;
payload_mass = IN.structures.payload_mass + IN.structures.extra_payload_mass;
payload_cg = 2;
tank_r = TANKS.singular.radius;
copv_mass = IN.press.copv.mass + 6;
copv_r = 0.18;
copv_t = 13.21/1000;
copv_h = 0.56;
copv_zcg = 1.21;
copv_ycg = 0.23;
it_h = TANKS.singular.radius;

[cg, MoI, ~, ~, ~, ~] = cg_moi_init(ox_mass, fu_mass, engine_mass, engine_cg, ox_tank_mass, ox_tank_cg, ox_tank_wall_thick, fu_tank_mass, fu_tank_cg, fu_tank_wall_thick, structures_mass, structures_cg, avi_mass, avi_cg, fluids_mass, fluids_cg, payload_mass, payload_cg, ox_tank_h, fu_tank_h, tank_r, it_h, copv_mass, copv_r, copv_t, copv_h, copv_zcg, copv_ycg);

ox_masses = linspace(0,IN.propulsion.oxidizer_mass,100);
fu_masses =  linspace(0,IN.propulsion.fuel_mass,100);

dcg_dMox_table = zeros(length(ox_masses), length(fu_masses));
dcg_dMfu_table = zeros(length(ox_masses), length(fu_masses));
dI_dMox_table = zeros(length(ox_masses), length(fu_masses), 9);
dI_dMfu_table = zeros(length(ox_masses), length(fu_masses), 9);

for i = 1:length(ox_masses)
    for j = 1:length(fu_masses)
        [cg, MoI, dcg_dMfu, dcg_dMox, dI_dMox, dI_dMfu] = cg_moi_init(ox_masses(i), fu_masses(j), engine_mass, engine_cg, ox_tank_mass, ox_tank_cg, ox_tank_wall_thick, fu_tank_mass, fu_tank_cg, fu_tank_wall_thick, structures_mass, structures_cg, avi_mass, avi_cg, fluids_mass, fluids_cg, payload_mass, payload_cg, ox_tank_h, fu_tank_h, tank_r, it_h, copv_mass, copv_r, copv_t, copv_h, copv_zcg, copv_ycg);
        dcg_dMox_table(i,j) = dcg_dMox;
        dcg_dMfu_table(i,j) = dcg_dMfu;
        dI_dMox_table(i,j, :) = reshape(dI_dMox,1,1,9);
        dI_dMfu_table(i,j, :) = reshape(dI_dMfu,1,1,9);
    end
end

save('cg_I_LUT.mat', ...
    "ox_masses", "fu_masses", ...
    "dcg_dMox_table", "dcg_dMfu_table", ...
    "dI_dMox_table", "dI_dMfu_table");
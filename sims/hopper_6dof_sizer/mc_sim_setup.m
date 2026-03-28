function mc_sim_setup(scenario)

% ---- Sizing (replaces main.m) ----
IN = mission_inputs(scenario);

TANKS  = size_tanks(IN);
AVI    = size_avionics(IN);
ENGINE = size_engine(IN);

VEH.mass.tanks   = TANKS.singular.total_mass;
VEH.mass.ox_tank = TANKS.singular.oxidizer.mass;
VEH.mass.fu_tank = TANKS.singular.fuel.mass;
VEH.h.ox_tank    = TANKS.singular.oxidizer.h;
VEH.h.fu_tank    = TANKS.singular.fuel.h;
VEH.mass.engine  = ENGINE.mass;
VEH.mass.avi     = AVI.mass;
VEH.mass.press   = IN.press.copv.mass * IN.press.copv.amount + IN.press.copv.gas_mass;
VEH.mass.dry     = IN.mass_factor * VEH.mass.engine + VEH.mass.tanks + VEH.mass.avi + ...
                   IN.structures.payload_mass + IN.structures.extra_payload_mass + VEH.mass.press;
VEH.mass.wet     = VEH.mass.dry + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;

[STRUCT, VEH] = size_structures(IN, VEH);

OUT = Outputs(IN, VEH, TANKS, STRUCT);

r_leg_1 = [0  IN.legs.feet_radial_dist  0];
r_leg_2 = [0 -IN.legs.feet_radial_dist  0];
r_leg_3 = [0  0  IN.legs.feet_radial_dist];
r_leg_4 = [0  0 -IN.legs.feet_radial_dist];

% ---- Propulsion ----
thrust_target = IN.propulsion.nominal_thrust;
ox_name       = IN.propulsion.oxidizer;
fu_name       = IN.propulsion.fuel;
fu_temp       = IN.propulsion.fu_temp;
ox_temp       = IN.propulsion.ox_temp;
ox_tank_P     = IN.propulsion.ox_press;
fu_tank_P     = IN.propulsion.fu_press;
ox_sys_CdA    = IN.propulsion.ox_sys_CdA;
fu_sys_CdA    = IN.propulsion.fu_sys_CdA;
eta_cstar     = IN.propulsion.eta_cstar;
At            = IN.propulsion.At;
eps           = IN.propulsion.eps;

[ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR, ox_valve_CdA, fu_valve_CdA, Isp] = ...
    prop_system(thrust_target, ox_name, fu_name, fu_temp, ox_temp, ...
                ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps);

% ---- Lookup tables (replaces load_lookup.m) ----
data = readmatrix('mdot_lookup.xlsx');
bp   = data(:,1);
tbl1 = data(:,2);
tbl2 = data(:,3);

load('wind_vectors2.mat');
load('cg_I_LUT.mat');

% ---- CG / MOI init (replaces cg_moi_test.m) ----
% Fixed geometry constants (from cg_moi_test.m)
engine_cg          = 0.3556;
ox_tank_cg         = 1.854;
ox_tank_wall_thick = 0.003175;
fu_tank_cg         = 1.0668;
fu_tank_wall_thick = 0.003175;
structures_cg      = 1.193;
avi_cg             = 1.4224;
fluids_cg          = 1.0668;
payload_cg         = 2;
copv_r             = 0.18;
copv_t             = 13.21/1000;
copv_h             = 0.56;
copv_zcg           = 1.21;
copv_ycg           = 0.23;

% Scenario/sized values
ox_mass        = IN.propulsion.oxidizer_mass;
fu_mass        = IN.propulsion.fuel_mass;
engine_mass    = ENGINE.mass;
ox_tank_mass   = TANKS.singular.oxidizer.mass;
ox_tank_h      = TANKS.singular.oxidizer.h;
fu_tank_mass   = TANKS.singular.fuel.mass;
fu_tank_h      = TANKS.singular.fuel.h;
structures_mass = OUT.Vehicle.MassDistribution.Structures;
avi_mass       = AVI.mass;
fluids_mass    = OUT.VehicleFluids.FinalMass;
payload_mass   = IN.structures.payload_mass + IN.structures.extra_payload_mass;
tank_r         = TANKS.singular.radius;
copv_mass      = IN.press.copv.mass + 6;
it_h           = TANKS.singular.radius;

[cg_init, MoI_init, ~, ~, ~, ~] = cg_moi_init( ...
    ox_mass, fu_mass, engine_mass, engine_cg, ...
    ox_tank_mass, ox_tank_cg, ox_tank_wall_thick, ...
    fu_tank_mass, fu_tank_cg, fu_tank_wall_thick, ...
    structures_mass, structures_cg, avi_mass, avi_cg, ...
    fluids_mass, fluids_cg, payload_mass, payload_cg, ...
    ox_tank_h, fu_tank_h, tank_r, it_h, ...
    copv_mass, copv_r, copv_t, copv_h, copv_zcg, copv_ycg);

% ---- MC scenario params (flat, for Simulink blocks) ----
slosh_lateral_damping   = scenario.slosh_lateral_damping;
slosh_axial_damping     = scenario.slosh_axial_damping;
throttle_rate_limit     = scenario.throttle_rate_limit * 4.44822;
throttle_latency        = scenario.throttle_latency;
tvc_actuator_rate_limit = deg2rad(scenario.tvc_actuator_rate_limit);
tvc_actuator_latency    = scenario.tvc_actuator_latency;
engine_off_axis_y       = scenario.engine_off_axis_y;
engine_off_axis_z       = scenario.engine_off_axis_z;
uwind = scenario.uwind;
vwind = scenario.vwind;
mode = 1
%azimuth                 = scenario.azimuth;

% ---- Push everything to base workspace (Simulink reads from here) ----
vars = who();
for k = 1:length(vars)
    assignin('base', vars{k}, eval(vars{k}));
end

end
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
% VEH.mass.dry     = IN.mass_factor * VEH.mass.engine + VEH.mass.tanks + VEH.mass.avi + ...
%                    IN.structures.payload_mass + IN.structures.extra_payload_mass + VEH.mass.press;

VEH.mass.dry = IN.mass_factor * 112.36;
VEH.mass.wet     = VEH.mass.dry + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;

[STRUCT, VEH] = size_structures(IN, VEH);

VEH.mass.dry = IN.mass_factor * 112.36;

VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass;

OUT = Outputs(IN, VEH, TANKS, STRUCT);

OUT.Vehicle.WetMass

VEH.mass.wet

r_leg_1 = [0  IN.legs.feet_radial_dist  0];
r_leg_2 = [0 -IN.legs.feet_radial_dist  0];
r_leg_3 = [0  0  IN.legs.feet_radial_dist];
r_leg_4 = [0  0 -IN.legs.feet_radial_dist];

VEH.k_engine = (IN.engine_cg/1.62)^2;

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
engine_mass    = 37.83;
ox_tank_mass   = TANKS.singular.oxidizer.mass;
ox_tank_h      = TANKS.singular.oxidizer.h;
fu_tank_mass   = TANKS.singular.fuel.mass;
fu_tank_h      = TANKS.singular.fuel.h;
structures_mass = 34.26;
avi_mass       = 11.52;
fluids_mass    = 28.41;
payload_mass   = IN.structures.payload_mass;
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
slosh_lateral_damping   = IN.tanks.lateral_damping;
slosh_axial_damping     = IN.tanks.axial_damping;
throttle_rate_limit     = IN.throttle_valve.rate_limit;
throttle_latency        = IN.throttle_valve_latency;
tvc_actuator_rate_limit = IN.tvc.max_gimbal_rate;
tvc_actuator_latency    = IN.tvc.actuator_latency;
engine_off_axis_y       = IN.off_axis_y;
engine_off_axis_z       = IN.off_axis_x;
uwind                   = IN.wind.uwind;
vwind                   = IN.wind.vwind;
wind_mode               = IN.wind.mode;
ox_mass                 = IN.propulsion.oxidizer_mass;
fu_mass                 = IN.propulsion.fuel_mass;

% ---- Push to base workspace ----
assignin('base', 'IN',                    IN);
assignin('base', 'VEH',                   VEH);
assignin('base', 'TANKS',                 TANKS);
assignin('base', 'STRUCT',                STRUCT);
assignin('base', 'OUT',                   OUT);
assignin('base', 'ENGINE',                ENGINE);
assignin('base', 'AVI',                   AVI);
assignin('base', 'r_leg_1',               r_leg_1);
assignin('base', 'r_leg_2',               r_leg_2);
assignin('base', 'r_leg_3',               r_leg_3);
assignin('base', 'r_leg_4',               r_leg_4);
assignin('base', 'ox_mdot',               ox_mdot);
assignin('base', 'fu_mdot',               fu_mdot);
assignin('base', 'tot_mdot',              tot_mdot);
assignin('base', 'Pc',                    Pc);
assignin('base', 'thrust',                thrust);
assignin('base', 'MR',                    MR);
assignin('base', 'ox_valve_CdA',          ox_valve_CdA);
assignin('base', 'fu_valve_CdA',          fu_valve_CdA);
assignin('base', 'Isp',                   Isp);
assignin('base', 'bp',                    bp);
assignin('base', 'tbl1',                  tbl1);
assignin('base', 'tbl2',                  tbl2);
assignin('base', 'engine_cg',             engine_cg);
assignin('base', 'cg_init',               cg_init);
assignin('base', 'MoI_init',              MoI_init);
assignin('base', 'ox_mass',               ox_mass);
assignin('base', 'fu_mass',               fu_mass);
assignin('base', 'tank_r',                tank_r);
assignin('base', 'slosh_lateral_damping', slosh_lateral_damping);
assignin('base', 'slosh_axial_damping',   slosh_axial_damping);
assignin('base', 'throttle_rate_limit',   throttle_rate_limit);
assignin('base', 'throttle_latency',      throttle_latency);
assignin('base', 'tvc_actuator_rate_limit', tvc_actuator_rate_limit);
assignin('base', 'tvc_actuator_latency',  tvc_actuator_latency);
assignin('base', 'engine_off_axis_y',     engine_off_axis_y);
assignin('base', 'engine_off_axis_z',     engine_off_axis_z);
assignin('base', 'uwind',                 uwind);
assignin('base', 'vwind',                 vwind);
assignin('base', 'wind_mode',             wind_mode);

end


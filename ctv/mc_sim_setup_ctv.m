function mc_sim_setup_ctv(scenario)
%MC_SIM_SETUP_CTV Build one CTV sim configuration in the base workspace.
%   mc_sim_setup_ctv(scenario) applies a Monte Carlo scenario struct.
%   mc_sim_setup_ctv([]) (or no argument) sets up the nominal case and
%   reproduces sim_setup_ctv.m exactly.
%
%   Same structure as sims/hopper_6dof_sizer/mc_sim_setup.m, but every
%   physics step is CTV's own: main_ctv sizing, load_lookup_ctv tables,
%   cg_moi_test_ctv CG/MoI, generateBaffleLUT_ctv slosh damping and
%   LinerizationMaster_ctv gains. The controller is designed on the nominal
%   dry block and the plant flies the dispersed one, so the mass dispersion
%   is something the controller has to reject.

persistent printed_which

if nargin < 1; scenario = []; end

% --- Setup CTV-only paths ---
% Deliberately NOT addpath(genpath(repoRoot)): sims/hopper_6dof_sizer holds
% same-named files (mission_inputs, Outputs, lateralDampingCalc, cg_moi_init,
% lookupMat, Trajectory2, ...) that would shadow CTV's.
ctvDir = fileparts(mfilename('fullpath'));
addpath(ctvDir);
addpath(fullfile(ctvDir, 'sizing'));
addpath(fullfile(ctvDir, 'inputs'));
addpath(fullfile(ctvDir, 'propulsion'));
addpath(fullfile(ctvDir, 'dynamics'));
addpath(fullfile(ctvDir, 'loads'));
addpath(fullfile(ctvDir, 'cg'));
addpath(fullfile(ctvDir, 'gnc'));
addpath(fullfile(ctvDir, 'baffle'));

% Print the resolved path of the files most likely to be shadowed, once per
% MATLAB process (so a shadowing bug is obvious instead of silent).
if isempty(printed_which)
    fprintf('mc_sim_setup_ctv path check:\n');
    fprintf('  mission_inputs     -> %s\n', which('mission_inputs'));
    fprintf('  lateralDampingCalc -> %s\n', which('lateralDampingCalc'));
    fprintf('  cg_moi_init        -> %s\n', which('cg_moi_init'));
    fprintf('  Outputs            -> %s\n', which('Outputs'));
    printed_which = true;
end

% ---- Sizing (main_ctv.m, minus the summary printout) ----
IN     = mission_inputs(scenario);
TANKS  = size_tanks(IN);
AVI    = size_avionics(IN);
ENGINE = size_engine(IN);

% Dry mass: the CTV hardcode from main_ctv.m, dispersed by mass_factor.
% size_structures then adds the landing legs and intertank on top of it,
% exactly as main_ctv.m lets it, which is why the nominal case reproduces
% sim_setup_ctv.m. VEH.mass.wet is the pre-legs value, also as in main_ctv.m.
dry_mass_nominal = 116.86;
VEH.mass.dry = IN.mass_factor * dry_mass_nominal;
VEH.mass.wet = VEH.mass.dry + ...
    IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;

[STRUCT, VEH] = size_structures(IN, VEH);

r_leg_1 = [0 IN.legs.feet_radial_dist 0];
r_leg_2 = [0 -IN.legs.feet_radial_dist 0];
r_leg_3 = [0 0 IN.legs.feet_radial_dist];
r_leg_4 = [0 0 -IN.legs.feet_radial_dist];

VEH.k_engine = (IN.engine_cg/IN.mission.max_vehicle_height)^2;

OUT = Outputs(IN, VEH, TANKS, STRUCT);

% Nominal dry block for the controller design: same sizing with mass_factor
% (and the cg/inertia factors) forced to 1. At nominal this is identical to
% OUT above, so the swap below is a no-op.
IN_nom = IN;
IN_nom.mass_factor    = 1.0;
IN_nom.cg_factor      = 1.0;
IN_nom.inertia_factor = 1.0;
VEH_nom.mass.dry = dry_mass_nominal;
VEH_nom.mass.wet = VEH_nom.mass.dry + ...
    IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;
[STRUCT_nom, VEH_nom] = size_structures(IN_nom, VEH_nom);
VEH_nom.k_engine = VEH.k_engine;
OUT_ctrl = Outputs(IN_nom, VEH_nom, TANKS, STRUCT_nom);

% ---- Lookup Tables & Spreadsheets (load_lookup_ctv.m) ----
mdot_data    = readmatrix('mdot_lookup.xlsx');
thrust_data  = mdot_data(:,1);   % breakpoints
fu_mdot_data = mdot_data(:,2);
ox_mdot_data = mdot_data(:,3);

% wind_vectors2.mat (uwndTable, vwndTable, pressureLevels) and cg_I_LUT.mat
% (ox_masses, fu_masses, dcg_dM*_table, dI_dM*_table): push everything they
% hold, so a new table in either file reaches the model without a code change.
for lutFile = {'wind_vectors2.mat', 'cg_I_LUT.mat'}
    lutData = load(lutFile{1});
    lutNames = fieldnames(lutData);
    for i = 1:length(lutNames)
        assignin('base', lutNames{i}, lutData.(lutNames{i}));
    end
end

% ---- CG / MOI Init (cg_moi_test_ctv.m) ----
ox_mass = IN.propulsion.oxidizer_mass;
fu_mass = IN.propulsion.fuel_mass;

engine_mass        = IN.engine_mass;
engine_cg          = IN.engine_cg;
mount_cg           = IN.mount_cg;
ox_tank_mass       = TANKS.singular.oxidizer.mass;
ox_tank_cg         = IN.ox_tank_cg;
ox_tank_wall_thick = IN.ox_tank_wall_thick;
ox_tank_h          = TANKS.singular.oxidizer.h;
fu_tank_mass       = TANKS.singular.fuel.mass;
fu_tank_cg         = IN.fu_tank_cg;
fu_tank_wall_thick = IN.fu_tank_wall_thick;
fu_tank_h          = TANKS.singular.fuel.h;
structures_mass    = IN.structures.structures_mass;
structures_cg      = IN.structures_cg;
avi_mass           = IN.avionics.avi_mass;
avi_cg             = IN.avi_cg;
fluids_mass        = IN.fluids_mass;
fluids_cg          = IN.fluids_cg;
payload_mass       = IN.structures.payload_mass;
payload_cg         = IN.payload_cg;
tank_r             = TANKS.singular.radius;
copv_mass          = IN.press.copv.mass + IN.press.copv.gas_mass;
copv_r             = IN.press.copv.r;
copv_t             = IN.press.copv.t;
copv_h             = IN.press.copv.h;
copv_zcg           = IN.copv_zcg;
copv_ycg           = IN.copv_ycg;
it_h               = TANKS.singular.radius;

[cg_init, MoI_init, ~, ~, ~, ~] = cg_moi_init(...
    ox_mass, fu_mass,...
    engine_mass, engine_cg, ...
    ox_tank_mass, ox_tank_cg, ox_tank_wall_thick, ...
    fu_tank_mass, fu_tank_cg, fu_tank_wall_thick, ...
    structures_mass, structures_cg, ...
    avi_mass, avi_cg,...
    fluids_mass, fluids_cg, ...
    payload_mass, payload_cg, ...
    ox_tank_h, fu_tank_h, tank_r, it_h, ...
    copv_mass, copv_r, copv_t, copv_h, copv_zcg, copv_ycg);

% CTV's cg_moi_init builds the vehicle from the component list and never sees
% the dry-mass hardcode, so mass_factor cannot move cg_init or MoI_init, and
% IN.cg_factor / IN.inertia_factor have nowhere to apply here (see
% mc_params_ctv.json, where they are deliberately absent). The controller and
% plant therefore share these two; the swap below is kept in main's shape so a
% real dry/propellant split can be dropped in later.
cg_init_ctrl  = cg_init;
MoI_init_ctrl = MoI_init;

% ---- Slosh damping tables (generateBaffleLUT_ctv.m) ----
% CTV's 7-argument damping model, evaluated at this scenario's propellant
% loads, then scaled per tank by the dispersed damping factor (log-uniform
% 0.5x-2x): the dispersion represents damping-model uncertainty, which dwarfs
% the effect of baffle manufacturing tolerances.
ox_density = 1141;   % kg/m^3
fu_density = 800;    % kg/m^3

fu_nu = 1.3e-6;  % m^2/s Fu
ox_nu = 0.17e-6; % m^2/s OX

ox_baffle_width  = TANKS.singular.radius * 0.3;
fu_baffle_width  = TANKS.singular.radius * 0.3;
ox_baffle_number = 5;
fu_baffle_number = 10;

ox_mass_profile = linspace(0.1, ox_mass, 200);
fu_mass_profile = linspace(0.1, fu_mass, 200);

ox_lateral_damping_ratios = IN.tanks.damping_factor_ox * lateralDampingCalc( ...
    ox_mass_profile, ox_density, ox_nu, tank_r, ox_baffle_width, ox_baffle_number, ox_tank_h);
fu_lateral_damping_ratios = IN.tanks.damping_factor_fu * lateralDampingCalc( ...
    fu_mass_profile, fu_density, fu_nu, tank_r, fu_baffle_width, fu_baffle_number, fu_tank_h);
ox_axial_damping_ratios   = IN.tanks.damping_factor_ox * axialDampingCalc( ...
    ox_mass_profile, ox_density, ox_nu, tank_r, ox_baffle_width, ox_baffle_number, ox_tank_h);
fu_axial_damping_ratios   = IN.tanks.damping_factor_fu * axialDampingCalc( ...
    fu_mass_profile, fu_density, fu_nu, tank_r, fu_baffle_width, fu_baffle_number, fu_tank_h);

% ---- MC Scenario Params (Flat for Simulink blocks) ----
% v2_ctv reads all of these through the IN struct rather than by these flat
% names, but they are pushed for parity with the 6DOF harness and so an MC run
% can be inspected from the base workspace. IN.tanks.lateral_damping and
% IN.tanks.axial_damping stay commented out in CTV, so unlike the 6DOF setup
% there is no slosh_lateral_damping / slosh_axial_damping here.
throttle_rate_limit     = IN.throttle_valve.rate_limit;
throttle_latency        = IN.throttle_valve_latency;
tvc_actuator_rate_limit = IN.tvc.max_gimbal_rate;
tvc_actuator_latency    = IN.tvc.actuator_latency;
engine_off_axis_y       = IN.off_axis_y;
engine_off_axis_z       = IN.off_axis_x;
uwind                   = IN.wind.uwind;
vwind                   = IN.wind.vwind;
wind_mode               = IN.wind.mode;

% ---- Run Linearization & Trajectory Generator ----
% LinerizationMaster_ctv and Trajectory2 are scripts that read and write the
% base workspace, so they are run there. They produce tgrid, unom, K1flat,
% K2grid, Target_Trajectory, the reference timeseries, g, and the Signal
% Editor dataset 'scenario' (which is why no MC scenario struct is ever
% assigned to a base variable of that name).
% The controller is designed on the NOMINAL dry block.
assignin('base', 'IN',        IN);
assignin('base', 'OUT',       OUT_ctrl);
assignin('base', 'VEH',       VEH_nom);
assignin('base', 'TANKS',     TANKS);
assignin('base', 'cg_init',   cg_init_ctrl);
assignin('base', 'MoI_init',  MoI_init_ctrl);
assignin('base', 'engine_cg', engine_cg);
evalin('base', 'LinerizationMaster_ctv');

% ---- Push Everything to Base Workspace for Simulink ----
% From here on the base workspace holds the DISPERSED plant.
assignin('base', 'IN',                        IN);
assignin('base', 'VEH',                       VEH);
assignin('base', 'TANKS',                     TANKS);
assignin('base', 'STRUCT',                    STRUCT);
assignin('base', 'OUT',                       OUT);
assignin('base', 'ENGINE',                    ENGINE);
assignin('base', 'AVI',                       AVI);
assignin('base', 'r_leg_1',                   r_leg_1);
assignin('base', 'r_leg_2',                   r_leg_2);
assignin('base', 'r_leg_3',                   r_leg_3);
assignin('base', 'r_leg_4',                   r_leg_4);
assignin('base', 'thrust_data',               thrust_data);
assignin('base', 'fu_mdot_data',              fu_mdot_data);
assignin('base', 'ox_mdot_data',              ox_mdot_data);
assignin('base', 'engine_cg',                 engine_cg);
assignin('base', 'engine_mass',               engine_mass);
assignin('base', 'mount_cg',                  mount_cg);
assignin('base', 'cg_init',                   cg_init);
assignin('base', 'MoI_init',                  MoI_init);
assignin('base', 'ox_mass',                   ox_mass);
assignin('base', 'fu_mass',                   fu_mass);
assignin('base', 'tank_r',                    tank_r);
assignin('base', 'ox_mass_profile',           ox_mass_profile);
assignin('base', 'fu_mass_profile',           fu_mass_profile);
assignin('base', 'ox_lateral_damping_ratios', ox_lateral_damping_ratios);
assignin('base', 'fu_lateral_damping_ratios', fu_lateral_damping_ratios);
assignin('base', 'ox_axial_damping_ratios',   ox_axial_damping_ratios);
assignin('base', 'fu_axial_damping_ratios',   fu_axial_damping_ratios);
assignin('base', 'throttle_rate_limit',       throttle_rate_limit);
assignin('base', 'throttle_latency',          throttle_latency);
assignin('base', 'tvc_actuator_rate_limit',   tvc_actuator_rate_limit);
assignin('base', 'tvc_actuator_latency',      tvc_actuator_latency);
assignin('base', 'engine_off_axis_y',         engine_off_axis_y);
assignin('base', 'engine_off_axis_z',         engine_off_axis_z);
assignin('base', 'uwind',                     uwind);
assignin('base', 'vwind',                     vwind);
assignin('base', 'wind_mode',                 wind_mode);

% ---- Verify the model's inputs are all there ----
% Every base-workspace variable Simulink.findVars reports for
% hopper_6dof_NED_v2_ctv. The six names findVars also reports but which live
% in block masks rather than the base workspace (const, flatIdxVect,
% flatLenVect, idxErrMode, numOutports, rowsOrCols) are excluded.
% 'g' is set to 9.81 by Trajectory2, not pushed here, so that the nominal case
% stays bit-identical to sim_setup_ctv.m.
required = { ...
    'ENGINE', 'IN', 'K1flat', 'K2grid', 'MoI_init', 'OUT', 'TANKS', 'T_ts', ...
    'Target_Trajectory', 'VEH', 'cg_init', 'dI_dMfu_table', 'dI_dMox_table', ...
    'dcg_dMfu_table', 'dcg_dMox_table', 'engine_cg', 'engine_mass', ...
    'fu_axial_damping_ratios', 'fu_lateral_damping_ratios', 'fu_mass_profile', ...
    'fu_masses', 'fu_mdot_data', 'g', 'mount_cg', 'ox_axial_damping_ratios', ...
    'ox_lateral_damping_ratios', 'ox_mass_profile', 'ox_masses', ...
    'ox_mdot_data', 'pressureLevels', 'tank_r', 'tgrid', 'thrust_data', ...
    'unom', 'uwndTable', 'vwndTable', 'vz_ref_ts', 'z_ref_ts', ...
    'STRUCT', 'AVI', 'throttle_rate_limit', 'throttle_latency', ...
    'tvc_actuator_rate_limit', 'tvc_actuator_latency', 'engine_off_axis_y', ...
    'engine_off_axis_z', 'uwind', 'vwind', 'wind_mode'};

inBase  = evalin('base', 'who');
missing = required(~ismember(required, inBase));
if ~isempty(missing)
    error('mc_sim_setup_ctv:missingBaseVars', ...
        'Setup finished but the base workspace is missing: %s', ...
        strjoin(missing, ', '));
end

% ---- Pin the random stream for the sim ----
% The model draws from MATLAB's global random stream (the wind model, which
% runs in Monte Carlo mode whenever a scenario is applied). Without this, the
% stream carries over between runs: repeated sims of the SAME scenario in one
% MATLAB process drifted by several m/s in landing velocity, and no amount of
% bdclose / clear functions / clear mex reset it. Seeding here makes each
% scenario reproducible and gives it its own wind realisation.
% generateScenarios_ctv already emits Seed = masterSeed + i for this.
% A nominal setup (empty scenario) deliberately leaves the stream untouched,
% so mc_sim_setup_ctv([]) still reproduces sim_setup_ctv.m exactly.
if ~isempty(scenario)
    if isfield(scenario, 'Seed')
        rng(scenario.Seed, 'twister');
    else
        rng(42, 'twister');   % master seed, for a scenario built from the JSON nominals
    end
end

end

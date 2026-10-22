function mc_sim_setup(scenario)
    % --- Setup Global Paths ---
    currentDir = fileparts(mfilename('fullpath'));
    addpath(currentDir);
    addpath(fullfile(currentDir, 'sizing'));
    addpath(fullfile(currentDir, 'inputs'));
    addpath(fullfile(currentDir, 'propulsion'));
    addpath(fullfile(currentDir, 'dynamics'));

    % ---- Sizing ----
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
    % Dry mass: main.m override (legs/intertank included), dispersed by mass_factor
    dry_mass_nominal = 116.86;
    VEH.mass.dry     = IN.mass_factor * dry_mass_nominal;
    VEH.mass.wet     = VEH.mass.dry + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;
    
    [STRUCT, VEH] = size_structures(IN, VEH);
    VEH.mass.dry = IN.mass_factor * dry_mass_nominal;   % size_structures adds legs; already included
    VEH.mass.wet = VEH.mass.dry + IN.propulsion.oxidizer_mass  + IN.propulsion.fuel_mass;
    
    OUT = Outputs(IN, VEH, TANKS, STRUCT);

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

    % ---- Lookup Tables & Spreadsheets ----
    mdot_data = readmatrix('mdot_lookup.xlsx');
    bp           = mdot_data(:,1);
    thrust_data  = mdot_data(:,1);   % Mapped for Simulink lookup blocks
    fu_mdot_data = mdot_data(:,2);
    ox_mdot_data = mdot_data(:,3);
    tbl1         = mdot_data(:,2);
    tbl2         = mdot_data(:,3);

    % ---- Dynamic .mat File Loading (Pushes ALL internal variables automatically) ----
    wind_data = load('wind_vectors2.mat');
    wind_fields = fieldnames(wind_data);
    for i = 1:length(wind_fields)
        assignin('base', wind_fields{i}, wind_data.(wind_fields{i}));
    end

    cg_data = load('cg_I_LUT.mat');
    cg_fields = fieldnames(cg_data);
    for i = 1:length(cg_fields)
        assignin('base', cg_fields{i}, cg_data.(cg_fields{i}));
    end

    % ---- CG / MOI Init ----
    engine_cg          = 0.3556;
    ox_tank_cg         = 1.854;
    ox_tank_wall_thick = 0.003175;
    fu_tank_cg         = 1.0668;
    fu_tank_wall_thick = 0.003175;
    structures_cg      = 1.193;
    avi_cg             = 1.4224;
    fluids_cg          = 1.0668;
    payload_cg         = 2;
    mount_cg           = IN.mount_cg;   % gimbal pivot (was the engine CG, 0.3556)
    copv_r             = 0.18;
    copv_t             = 13.21/1000;
    copv_h             = 0.56;
    copv_zcg           = 1.21;
    copv_ycg           = 0.23;

    ox_mass         = IN.propulsion.oxidizer_mass;
    fu_mass         = IN.propulsion.fuel_mass;
    engine_mass     = 37.83;
    ox_tank_mass    = TANKS.singular.oxidizer.mass;
    ox_tank_h       = IN.tanks.ox_h_as_built;   % as-built, not sizer TANKS.singular.oxidizer.h
    fu_tank_mass    = TANKS.singular.fuel.mass;
    fu_tank_h       = IN.tanks.fu_h_as_built;   % as-built, not sizer TANKS.singular.fuel.h
    structures_mass = 34.26;
    avi_mass        = 11.52;
    fluids_mass     = 28.41;
    payload_mass    = IN.structures.payload_mass;
    tank_r          = TANKS.singular.radius;
    copv_mass       = IN.press.copv.mass + 6;
    it_h            = TANKS.singular.radius;

    % Dry block plus propellant columns. The plant gets this scenario's
    % dispersed dry mass, CG and inertia; the controller is designed on the
    % nominal dry block (swapped in around LinerizationMaster below), so the
    % dispersion is something the controller has to reject.
    MP = mass_properties_setup(IN, TANKS, VEH.mass.dry);
    IN_nom = IN; IN_nom.cg_factor = 1; IN_nom.inertia_factor = 1;
    MP_ctrl = mass_properties_setup(IN_nom, TANKS, dry_mass_nominal);
    OUT_ctrl = OUT;
    OUT_ctrl.Vehicle.DryMass = dry_mass_nominal;
    OUT_ctrl.Vehicle.WetMass = dry_mass_nominal + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;
    cg_init  = MP.cg_init;
    MoI_init = MP.MoI_init;
    for f = {'ox_masses', 'fu_masses', 'dcg_dMox_table', 'dcg_dMfu_table', 'dI_dMox_table', 'dI_dMfu_table'}
        assignin('base', f{1}, MP.(f{1}));
    end

    % ---- Slosh & Profile Dummies (Safeguards for Lookup Blocks) ----
    ox_mass_profile         = linspace(0, ox_mass, 100); % Fallback profile if needed
    ox_lateral_damping_ratios = IN.tanks.lateral_damping * ones(size(ox_mass_profile));
    ox_axial_damping_ratios   = IN.tanks.axial_damping   * ones(size(ox_mass_profile));

    fu_mass_profile           = linspace(0, fu_mass, 100);                    
    fu_lateral_damping_ratios = IN.tanks.lateral_damping * ones(size(fu_mass_profile)); 
    fu_axial_damping_ratios   = IN.tanks.axial_damping   * ones(size(fu_mass_profile));

    % ---- MC Scenario Params (Flat for Simulink blocks) ----
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

    % ---- Run Linearization & Trajectory Generator ----
    % This generates tgrid, unom, K1flat, K2grid, Target_Trajectory, reference timeseries, etc.
    OUT_plant = OUT;     MP_plant = MP;
    OUT = OUT_ctrl;      MP = MP_ctrl;
    cg_init = MP.cg_init;  MoI_init = MP.MoI_init;
    LinerizationMaster();
    OUT = OUT_plant;     MP = MP_plant;
    cg_init = MP.cg_init;  MoI_init = MP.MoI_init;

    % ---- Push Everything to Base Workspace for Simulink ----
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
    assignin('base', 'thrust_data',           thrust_data);
    assignin('base', 'fu_mdot_data',          fu_mdot_data);
    assignin('base', 'ox_mdot_data',          ox_mdot_data);
    assignin('base', 'ox_mass_profile',       ox_mass_profile);
    assignin('base', 'ox_lateral_damping_ratios', ox_lateral_damping_ratios);
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
    assignin('base', 'engine_mass',           engine_mass);
    assignin('base', 'g',                     IN.const.g0);
    assignin('base', 'mount_cg', mount_cg);
    assignin('base', 'fu_mass_profile',           fu_mass_profile);         
    assignin('base', 'fu_lateral_damping_ratios', fu_lateral_damping_ratios);  
    assignin('base', 'ox_axial_damping_ratios',   ox_axial_damping_ratios);   
    assignin('base', 'fu_axial_damping_ratios',   fu_axial_damping_ratios);  
end
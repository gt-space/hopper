function results = mc_runner(n_scenarios, params_file, output_file)

if nargin < 1; n_scenarios = 100;              end
if nargin < 2; params_file = 'mc_params.json'; end
if nargin < 3; output_file = 'mc_results.mat'; end

addpath('./sizing')
addpath('./inputs')
addpath('./propulsion')
addpath('./dynamics')

load_system('hopper_6dof_NED_v2');
% set_param('hopper_6dof_NED_v2', 'SolverType', 'Fixed-step')
% set_param('hopper_6dof_NED_v2', 'FixedStep', '0.0005')

scenarios    = generateScenarios(params_file, n_scenarios);
assignin('base', 'scenarios',             scenarios);
n_scenarios  = length(scenarios);
results_cell = cell(n_scenarios, 1);
t_start      = tic;

fprintf('Starting Monte Carlo: %d scenarios\n', n_scenarios);

fprintf('Running nominal scenario...\n');
nominal = runNominal(params_file);


for mc_iter = 1:n_scenarios
    fprintf('Running scenario %d / %d\n', mc_iter, n_scenarios);
    try
        mc_sim_setup(scenarios(mc_iter));
        IN        = evalin('base', 'IN');
        VEH       = evalin('base', 'VEH');
        TANKS     = evalin('base', 'TANKS');
        STRUCT    = evalin('base', 'STRUCT');
        cg_init   = evalin('base', 'cg_init');
        engine_cg = evalin('base', 'engine_cg');
        OUT       = Outputs(IN, VEH, TANKS, STRUCT);
        LinerizationMaster();
        mc_sim_setup(scenarios(mc_iter));
        sim_out = sim('hopper_6dof_NED_v2');
        r = mc_main(scenarios(mc_iter), sim_out);
        r = check_constraints(r);
    catch err
        fprintf('Run %d error: %s\n', mc_iter, err.message);
        r = make_empty_result(scenarios(mc_iter), err.message);
    end
    results_cell{mc_iter} = r;
end

results = [results_cell{:}];

n_pass  = sum(arrayfun(@(r) isfield(r,'status') && r.status.success && r.status.pass,  results));
n_fail  = sum(arrayfun(@(r) isfield(r,'status') && r.status.success && ~r.status.pass, results));
n_error = sum(arrayfun(@(r) isfield(r,'status') && ~r.status.success,                  results));

fprintf('%-6s %-8s %-10s %-10s %-10s %-6s\n', ...
    'Run', 'Status', 'DryMass', 'WetMass', 'TWR_init', 'Pass');
fprintf('%s\n', repmat('-', 1, 52));

for i = 1:n_scenarios
    r = results(i);
    if ~r.status.success
        fprintf('%-6d %-8s  %s\n', i, 'ERROR', r.status.error);
    elseif r.status.pass
        fprintf('%-6d %-8s %-10.2f %-10.2f %-10.3f %-6s\n', ...
            i, 'OK', r.vehicle.dry_mass, r.vehicle.wet_mass, r.vehicle.twr_initial, 'PASS');
    else
        fprintf('%-6d %-8s %-10.2f %-10.2f %-10.3f %-6s\n', ...
            i, 'OK', r.vehicle.dry_mass, r.vehicle.wet_mass, r.vehicle.twr_initial, 'FAIL');
    end
end

elapsed = toc(t_start);

fprintf('\n%s\n', repmat('=', 1, 52));
fprintf('Complete:  %d scenarios in %.1f s\n', n_scenarios, elapsed);
fprintf('Pass:      %d (%.1f%%)\n', n_pass,  100*n_pass/n_scenarios);
fprintf('Fail:      %d (%.1f%%)\n', n_fail,  100*n_fail/n_scenarios);
fprintf('Errors:    %d (%.1f%%)\n', n_error, 100*n_error/n_scenarios);
fprintf('%s\n', repmat('=', 1, 52));

save(output_file, 'results', 'scenarios', 'nominal');  % line ~65
fprintf('Results saved to %s\n', output_file);

end

function r = make_empty_result(scenario, error_msg)
r.scenario         = scenario;
r.vehicle.dry_mass    = NaN;
r.vehicle.wet_mass    = NaN;
r.vehicle.mass_ratio  = NaN;
r.vehicle.twr_initial = NaN;
r.vehicle.twr_final   = NaN;
r.vehicle.thrust_min  = NaN;
r.vehicle.thrust_max  = NaN;
r.propellant.ox_mass      = NaN;
r.propellant.fuel_mass    = NaN;
r.propellant.ox_press_psi = NaN;
r.propellant.fu_press_psi = NaN;
r.tanks.ox_volume_L = NaN;
r.tanks.fu_volume_L = NaN;
r.tanks.ox_mass     = NaN;
r.tanks.ox_height   = NaN;
r.tanks.fu_mass     = NaN;
r.tanks.fu_height   = NaN;
r.tanks.total_mass  = NaN;
r.tanks.radius      = NaN;
r.copv.mass          = NaN;
r.copv.gas_mass      = NaN;
r.copv.amount        = NaN;
r.copv.volume        = NaN;
r.copv.max_press_psi = NaN;
r.structures.mass    = NaN;
r.avionics.battery_capacity = NaN;
r.avionics.num_cells        = NaN;
r.avionics.mass             = NaN;
r.engine.eta_cstar   = NaN;
r.engine.throat_area = NaN;
r.engine.eps         = NaN;
r.engine.inj_mass    = NaN;
r.engine.TCA_mass    = NaN;
r.engine.TVC_mass    = NaN;
r.engine.total_mass  = NaN;
r.flight.max_altitude    = NaN;
r.flight.max_ascent_vel  = NaN;
r.flight.max_descent_vel = NaN;
r.flight.landing_vel     = NaN;
r.flight.flight_time     = NaN;
r.flight.final_ox_mass   = NaN;
r.flight.final_fu_mass   = NaN;
r.flight.cg_init         = NaN;
r.flight.cg_final        = NaN;
r.timeseries.t        = [];
r.timeseries.x        = [];
r.timeseries.y        = [];
r.timeseries.z        = [];
r.timeseries.ox_mass  = [];
r.timeseries.fu_mass  = [];
r.timeseries.tot_mass = [];
r.timeseries.vz       = [];
r.status.success = false;
r.status.error   = error_msg;
r.status.pass    = false;
end

function result = check_constraints(result)

v = result.vehicle;
f = result.flight;
passes = true;

if v.twr_initial        < 1.0; passes = false; end
if v.twr_final          < 1.0; passes = false; end
if v.dry_mass           <= 0;  passes = false; end
if v.wet_mass           <= 0;  passes = false; end
if f.max_altitude       < 51;  passes = false; end
if f.max_ascent_vel     > 7;   passes = false; end
if f.max_descent_vel    > 7;   passes = false; end
if f.landing_vel        > 0.5; passes = false; end

result.status.pass = passes;



end

function nominal = runNominal(params_file)
    fid  = fopen(params_file);
    raw  = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    params = jsondecode(raw);
    fields = fieldnames(params);
    
    % Build nominal scenario from JSON nominal values
    scenario = struct();
    for k = 1:length(fields)
        scenario.(fields{k}) = params.(fields{k}).nominal;
    end
    
    mc_sim_setup(scenario);
    IN        = evalin('base', 'IN');
    VEH       = evalin('base', 'VEH');
    TANKS     = evalin('base', 'TANKS');
    STRUCT    = evalin('base', 'STRUCT');
    cg_init   = evalin('base', 'cg_init');
    engine_cg = evalin('base', 'engine_cg');
    OUT       = Outputs(IN, VEH, TANKS, STRUCT);
    saved_scenario = scenario;
    LinerizationMaster();
    scenario  = saved_scenario;
    mc_sim_setup(scenario);
    sim_out   = sim('hopper_6dof_NED_v2');
    nominal   = mc_main(scenario, sim_out);
end
function results = mc_runner(n_scenarios, params_file, output_file)

if nargin < 1; n_scenarios = 500;              end
if nargin < 2; params_file = 'mc_params.json'; end
if nargin < 3; output_file = 'mc_results.mat'; end

addpath('./sizing')
addpath('./inputs')
addpath('./propulsion')
addpath('./dynamics')

scenarios = generateScenarios(params_file, n_scenarios);
results_cell = cell(n_scenarios, 1);

n_pass  = 0;
n_fail  = 0;
n_error = 0;
t_start = tic;

fprintf('Starting Monte Carlo: %d scenarios\n', n_scenarios);

parfor i = 1:n_scenarios
    try
        r = mc_main(scenarios(i));
        r = check_constraints(r);
    catch err
        r                   = struct();
        r.scenario          = scenarios(i);
        r.status.success    = false;
        r.status.error      = err.message;
        r.status.pass       = false;
    end
    results_cell{i} = r;
end

results = [results_cell{:}];

fprintf('%-6s %-8s %-10s %-10s %-10s %-6s\n', ...
    'Run', 'Status', 'DryMass', 'WetMass', 'TWR_init', 'Pass');
fprintf('%s\n', repmat('-', 1, 52));

for i = 1:n_scenarios
    r = results(i);
    if ~r.status.success
        n_error = n_error + 1;
        fprintf('%-6d %-8s  %s\n', i, 'ERROR', r.status.error);
    elseif r.status.pass
        n_pass = n_pass + 1;
        fprintf('%-6d %-8s %-10.2f %-10.2f %-10.3f %-6s\n', ...
            i, 'OK', r.vehicle.dry_mass, r.vehicle.wet_mass, r.vehicle.twr_initial, 'PASS');
    else
        n_fail = n_fail + 1;
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

save(output_file, 'results', 'scenarios');
fprintf('Results saved to %s\n', output_file);

end


function result = check_constraints(result)

v = result.vehicle;
passes = true;

if v.twr_initial < 1.0; passes = false; end
if v.twr_final   < 1.0; passes = false; end
if v.dry_mass    <= 0;  passes = false; end
if v.wet_mass    <= 0;  passes = false; end

result.status.pass = passes;

end
% --- Monte Carlo Simulation Master Script ---

warning('off', 'MATLAB:Python:PyNotFound')
clear; clc;

n = input('Enter the number of Monte Carlo scenarios (e.g., 1000): ');
if isempty(n) || n <= 0
    error('Invalid input. Please enter a positive integer.');
end

if isempty(gcp('nocreate'))
    parpool();
end

% Attach all required models, lookup files, and subfolder directories to workers
p = gcp();
addAttachedFiles(p, { ...
    'hopper_6dof_NED_v2.slx', ...
    'mc_params.json', ...
    'mdot_lookup.xlsx', ...
    'wind_vectors2.mat', ...
    'cg_I_LUT.mat' ...
    });

jsonFile = 'mc_params.json';
[scenarios, mcTable] = generateScenarios(jsonFile, n);

scenarioStructs = table2struct(mcTable);
resultsCell = cell(n, 1); 

% --- Parallel Execution Loop ---
tic;
parfor i = 1:n
    % Ensure all helper directories are visible on this worker thread
    addpath(genpath(pwd)); 

    currentScenario = scenarioStructs(i);
    localResult = struct(); % Initialize local loop variable for classification

    try
        % Execute setup function to populate the worker's base workspace
        mc_sim_setup(currentScenario);

        % Configure Simulink simulation input object
        simInput = Simulink.SimulationInput('hopper_6dof_NED_v2');
        simInput = simInput.setVariable('currentScenario', currentScenario);

        % Execute simulation
        sim_out = sim(simInput);

        % Post-process results using mc_main
        localResult = mc_main(currentScenario, sim_out); 
    catch ME
        % Handle simulation crashes gracefully within the local scope
        localResult.scenario = currentScenario;
        localResult.status.success = false;
        localResult.status.error = ME.message;
    end

    resultsCell{i} = localResult;
end
elapsedTime = toc;
fprintf('Completed %d Monte Carlo runs in %.2f seconds.\n', n, elapsedTime);

% --- Export Results ---
mcTable.Results = resultsCell;
writetable(mcTable, 'mc_results_parallel.csv');
fprintf('Results successfully saved to mc_results_parallel.csv\n');
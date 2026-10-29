% --- Monte Carlo Simulation Master Script (Per-Run Checkpointed) ---
warning('off', 'MATLAB:Python:PyNotFound')
clear; clc;

% --- Setup Global Paths ---
currentDir = fileparts(mfilename('fullpath'));
addpath(genpath(currentDir)); % Recursively adds all subfolders automatically
addpath(fullfile(currentDir, 'sizing'));
addpath(fullfile(currentDir, 'inputs'));
addpath(fullfile(currentDir, 'propulsion'));
addpath(fullfile(currentDir, 'dynamics'));

% Options: 'parallel', 'serial', 'nominal'
runMode = 'parallel';

if strcmp(runMode, 'parallel')
    % ------------------------------------------------------------
    % User-tunable memory settings
    % ------------------------------------------------------------
    RAM_per_worker_GB = 4.0;   % Estimated RAM required per worker
    RAM_reserve_GB    = 6.0;   % RAM reserved for OS + MATLAB client
    
    % ------------------------------------------------------------
    % Get available system RAM
    % ------------------------------------------------------------
    try
        [~, sys] = memory;
        available_RAM_GB = sys.PhysicalMemory.Available / 1024^3;
    catch
        warning('Could not determine available system RAM.');
        warning('Using 1 worker as a safe fallback.');
        available_RAM_GB = RAM_reserve_GB + RAM_per_worker_GB;
    end
    
    % ------------------------------------------------------------
    % Calculate worker limits
    % ------------------------------------------------------------
    usable_RAM_GB = max(0, available_RAM_GB - RAM_reserve_GB);
    workers_by_RAM = floor(usable_RAM_GB / RAM_per_worker_GB);
    workers_by_CPU = feature('numcores');
    
    numWorkers = min([workers_by_RAM, workers_by_CPU]);
    numWorkers = max(1, numWorkers);
    
    fprintf('\n');
    fprintf('===============================================\n');
    fprintf(' Dynamic Parallel Worker Configuration\n');
    fprintf('===============================================\n');
    fprintf('Available RAM:       %.2f GB\n', available_RAM_GB);
    fprintf('Reserved RAM:        %.2f GB\n', RAM_reserve_GB);
    fprintf('Estimated RAM/worker: %.2f GB\n', RAM_per_worker_GB);
    fprintf('RAM worker limit:    %d\n', workers_by_RAM);
    fprintf('CPU worker limit:    %d\n', workers_by_CPU);
    fprintf('Selected workers:    %d\n', numWorkers);
    fprintf('===============================================\n\n');
    
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end
    parpool('local', numWorkers);
end

% Attach all required models, lookup files, and subfolder directories to workers
mFiles = dir('*.m');
matFiles = dir('*.mat');
slxFiles = dir('*.slx');
xlsxFiles = dir('*.xlsx');
jsonFiles = dir('*.json');
addAllFiles = [{matFiles.name}, {slxFiles.name}, {xlsxFiles.name}, {jsonFiles.name}, {mFiles.name}];
if strcmp(runMode, 'parallel')
    p = gcp();
    addAttachedFiles(p, addAllFiles);
end

jsonFile = 'mc_params.json';
ckptDir = fullfile(currentDir, 'mc_checkpoints'); % Use absolute path

if ~exist(ckptDir, 'dir')
    mkdir(ckptDir);
end

if ~strcmp(runMode, 'nominal')
    n = input('Enter the number of Monte Carlo scenarios (e.g., 1000): ');
    if isempty(n) || n <= 0
        error('Invalid input. Please enter a positive integer.');
    end
    
    [scenarios, mcTable] = generateScenarios(jsonFile, n);
    scenarioStructs = table2struct(mcTable);
    
    % --- Load existing checkpoints per scenario ---
    resultsCell = cell(n, 1);
    for i = 1:n
        ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', i));
        if exist(ckptFile, 'file')
            data = load(ckptFile, 'localResult');
            resultsCell{i} = data.localResult;
        end
    end
    completedCount = sum(~cellfun(@isempty, resultsCell));
    fprintf('Checkpoint status: %d / %d scenarios already completed.\n', completedCount, n);
end

% --- Execution Router (Switch-Case) ---
switch runMode
    case 'parallel'
        tic;
        missingIdx = find(cellfun(@isempty, resultsCell));
        numMissing = length(missingIdx);
        
        if numMissing == 0
            fprintf('All %d scenarios are already completed from checkpoints!\n', n);
        else
            fprintf('Preparing %d remaining scenarios for parallel execution...\n', numMissing);
            
            simIn = repmat(Simulink.SimulationInput('hopper_6dof_NED_v2'), numMissing, 1);
            for k = 1:numMissing
                i = missingIdx(k);
                currentScenario = scenarioStructs(i);
                
                simIn(k) = simIn(k).setVariable('currentScenario', currentScenario);
                simIn(k) = simIn(k).setPreSimFcn(@(in) local_pre_sim(in, currentScenario));
                
                % Pass absolute ckptDir into the post-sim function
                simIn(k) = simIn(k).setPostSimFcn(@(out) local_post_sim(out, currentScenario, i, ckptDir));
            end
            
            fprintf('Running parsim across workers...\n');
            simOuts = parsim(simIn, 'ShowProgress', 'on');
            
            % Reload the newly generated individual checkpoint files into resultsCell
            for k = 1:numMissing
                i = missingIdx(k);
                ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', i));
                if exist(ckptFile, 'file')
                    data = load(ckptFile, 'localResult');
                    resultsCell{i} = data.localResult;
                end
            end
        end
        
        elapsedTime = toc;
        fprintf('Completed parallel execution phase in %.2f seconds.\n', elapsedTime);
        
        % --- Run Nominal Case for Comparison ---
        fprintf('Running nominal case...\n');
        nominalScenario = mc_inputs();
        mc_sim_setup(nominalScenario);
        simInputNom = Simulink.SimulationInput('hopper_6dof_NED_v2');
        simOutNom = sim(simInputNom);
        nominal = mc_main(nominalScenario, simOutNom);
        
        % --- Export Final Results ---
        results = [resultsCell{:}];
        save('mc_results_parallel.mat', 'results', 'nominal');
        fprintf('Final results and nominal data successfully saved to mc_results_parallel.mat\n');
        
    case 'serial'
        tic;
        for i = 1:n
            if ~isempty(resultsCell{i})
                continue;
            end
            
            addpath(genpath(pwd)); 
            addpath('./sizing'); addpath('./inputs'); addpath('./propulsion'); addpath('./dynamics');
            
            currentScenario = scenarioStructs(i);
            localResult = struct(); 
            
            try
                mc_sim_setup(currentScenario);
                simInput = Simulink.SimulationInput('hopper_6dof_NED_v2');
                simInput = simInput.setVariable('currentScenario', currentScenario);
                
                sim_out = sim(simInput);
                localResult = mc_main(currentScenario, sim_out); 
                localResult = check_constraints(localResult);
            catch ME
                localResult.scenario = currentScenario;
                localResult.status.success = false;
                localResult.status.error = ME.message;
                localResult.status.pass = false;
            end
            
            resultsCell{i} = localResult;
            
            % Save individual checkpoint file immediately after each run
            ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', i));
            save(ckptFile, 'localResult');
        end
        
        elapsedTime = toc;
        fprintf('Completed sequential execution runs in %.2f seconds.\n', elapsedTime);
        
        mcTable.Results = resultsCell;
        save('mc_results_serial.mat', 'mcTable');
        fprintf('Final results successfully saved to mc_results_serial.mat\n');
        
   case 'nominal'
        tic;
        try
            addpath(genpath(pwd)); 
            addpath('./sizing'); addpath('./inputs'); addpath('./propulsion'); addpath('./dynamics');
            
            nominalScenario.ox_mass                 = 16;
            nominalScenario.fuel_mass               = 13;
            nominalScenario.cstar                   = 0.85;
            nominalScenario.mass_factor             = 1.0;
            nominalScenario.slosh_lateral_damping   = 0.08; 
            nominalScenario.slosh_axial_damping     = 0.3;  
            nominalScenario.throttle_rate_limit     = 180;  
            nominalScenario.throttle_latency        = 0.01;
            nominalScenario.tvc_actuator_rate_limit = 15; 
            nominalScenario.tvc_actuator_latency    = 0.01;
            nominalScenario.engine_off_axis_y       = 0;
            nominalScenario.engine_off_axis_z       = 0;
            nominalScenario.uwind                   = 5;
            nominalScenario.vwind                   = 5;
            
            mc_sim_setup(nominalScenario);
            
            simInput = Simulink.SimulationInput('hopper_6dof_NED_v2');
            sim_out = sim(simInput);
            nominal = mc_main(nominalScenario, sim_out);
            
            save('mc_results_nominal.mat', 'nominal');
            disp('Nominal simulation completed and saved successfully.');
            
        catch ME
            rethrow(ME);
        end
        elapsedTime = toc;
        fprintf('Nominal run completed in %.2f seconds.\n', elapsedTime);
        
    otherwise
        error('Invalid runMode specified. Use ''parallel'', ''serial'', or ''nominal''.');
end

function in = local_pre_sim(in, scenario)
    addpath(genpath(pwd));
    mc_sim_setup(scenario);
end

function out = local_post_sim(out, scenario, idx, ckptDir)
    % This runs on the worker thread the exact moment a simulation finishes
    addpath(genpath(pwd));
    try
        if ~isempty(out.ErrorMessage)
            error(out.ErrorMessage);
        end
        localResult = mc_main(scenario, out);
    catch ME
        localResult.scenario = scenario;
        localResult.status.success = false;
        localResult.status.error = ME.message;
        localResult.status.pass = false;
    end
    
    % Save individual checkpoint file instantly to the absolute path
    ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', idx));
    save(ckptFile, 'localResult');
end
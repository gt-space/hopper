% --- Monte Carlo Simulation Master Script (Auto-Batching & Checkpointed) ---
warning('off', 'MATLAB:Python:PyNotFound')
clear; clc;

% --- Setup Global Paths ---
currentDir = fileparts(mfilename('fullpath'));
addpath(genpath(currentDir)); % Recursively adds all subfolders automatically
addpath(fullfile(currentDir, 'sizing'));
addpath(fullfile(currentDir, 'inputs'));
addpath(fullfile(currentDir, 'propulsion'));
addpath(fullfile(currentDir, 'dynamics'));

% --- Execution & Batch Configuration ---
runMode     = 'parallel'; % Options: 'parallel', 'serial', 'nominal'
useBatching = true;       % Enable batch file splitting
batchSize   = 300;        % Number of scenarios per batch file (e.g., 300)

% --- Setup Directories ---
jsonFile  = 'mc_params.json';
ckptDir   = fullfile(currentDir, 'mc_checkpoints'); % Individual scenario checkpoints
batchDir  = fullfile(currentDir, 'mc_batches');     % Batch definition files

if ~exist(ckptDir, 'dir'), mkdir(ckptDir); end
if ~exist(batchDir, 'dir'), mkdir(batchDir); end

% --- Parallel Pool Setup ---
if strcmp(runMode, 'parallel')
    RAM_per_worker_GB = 4.0;   
    RAM_reserve_GB    = 6.0;   
    
    try
        [~, sys] = memory;
        available_RAM_GB = sys.PhysicalMemory.Available / 1024^3;
    catch
        warning('Could not determine available system RAM. Using 1 worker as fallback.');
        available_RAM_GB = RAM_reserve_GB + RAM_per_worker_GB;
    end
    
    usable_RAM_GB = max(0, available_RAM_GB - RAM_reserve_GB);
    workers_by_RAM = floor(usable_RAM_GB / RAM_per_worker_GB);
    workers_by_CPU = feature('numcores');
    
    numWorkers = min([workers_by_RAM, workers_by_CPU]);
    numWorkers = max(1, numWorkers);
    
    fprintf('\n===============================================\n');
    fprintf(' Dynamic Parallel Worker Configuration\n');
    fprintf('===============================================\n');
    fprintf('Available RAM:       %.2f GB\n', available_RAM_GB);
    fprintf('Selected workers:    %d\n', numWorkers);
    fprintf('===============================================\n\n');
    
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end
    parpool('local', numWorkers);
end

% Attach all required models and files to workers
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

% --- Scenario & Batch Generation / Loading ---
if ~strcmp(runMode, 'nominal')
    masterScenariosFile = fullfile(batchDir, 'master_scenarios.mat');
    
    if ~exist(masterScenariosFile, 'file')
        n = input('Enter the total number of Monte Carlo scenarios (e.g., 1000): ');
        if isempty(n) || n <= 0
            error('Invalid input. Please enter a positive integer.');
        end
        [scenarios, mcTable] = generateScenarios(jsonFile, n);
        
        % Split into batches
        numBatches = ceil(n / batchSize);
        allScenarioStructs = table2struct(mcTable);
        
        for b = 1:numBatches
            startIdx = (b - 1) * batchSize + 1;
            endIdx   = min(b * batchSize, n);
            
            batchData.batchID     = b;
            batchData.startIdx    = startIdx;
            batchData.endIdx      = endIdx;
            batchData.scenarios   = allScenarioStructs(startIdx:endIdx);
            batchData.mcTable     = mcTable(startIdx:endIdx, :);
            
            batchFile = fullfile(batchDir, sprintf('batch_%03d.mat', b));
            save(batchFile, 'batchData');
        end
        save(masterScenariosFile, 'n', 'numBatches');
        fprintf('Successfully generated and split %d scenarios into %d batches of size %d.\n', n, numBatches, batchSize);
    else
        load(masterScenariosFile, 'n', 'numBatches');
    end
    
    % --- AUTOMATICALLY FIND THE NEXT INCOMPLETE BATCH ---
    targetBatch = [];
    for b = 1:numBatches
        batchFile = fullfile(batchDir, sprintf('batch_%03d.mat', b));
        bData = load(batchFile, 'batchData');
        gIndices = bData.batchData.startIdx : bData.batchData.endIdx;
        
        % Check if all checkpoints exist for this batch
        batchComplete = true;
        for idx = gIndices
            ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', idx));
            if ~exist(ckptFile, 'file')
                batchComplete = false;
                break;
            end
        end
        
        if ~batchComplete
            targetBatch = b;
            break;
        end
    end
    
    if isempty(targetBatch)
        fprintf('\n===============================================\n');
        fprintf(' 🎉 All %d scenarios across all batches are fully completed!\n', n);
        fprintf('===============================================\n');
        return;
    end
    
    % --- Load Target Batch ---
    batchFile = fullfile(batchDir, sprintf('batch_%03d.mat', targetBatch));
    loadedBatch = load(batchFile, 'batchData');
    batchData   = loadedBatch.batchData;
    
    scenarioStructs = batchData.scenarios;
    globalIndices   = batchData.startIdx : batchData.endIdx;
    batchNumScen    = length(scenarioStructs);
    
    fprintf('\n[Auto-Selected] Running Batch %d/%d (Global Scenarios %d to %d).\n', targetBatch, numBatches, batchData.startIdx, batchData.endIdx);
    
    % --- Load existing checkpoints for this specific batch ---
    resultsCell = cell(batchNumScen, 1);
    for k = 1:batchNumScen
        i = globalIndices(k); 
        ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', i));
        if exist(ckptFile, 'file')
            data = load(ckptFile, 'localResult');
            resultsCell{k} = data.localResult;
        end
    end
    completedCount = sum(~cellfun(@isempty, resultsCell));
    fprintf('Batch checkpoint status: %d / %d scenarios already completed.\n', completedCount, batchNumScen);
end

% --- Execution Router (Switch-Case) ---
switch runMode
    case 'parallel'
        tic;
        missingLocalIdx = find(cellfun(@isempty, resultsCell));
        numMissing = length(missingLocalIdx);
        
        if numMissing == 0
            fprintf('All scenarios in Batch %d are already completed from checkpoints!\n', targetBatch);
        else
            fprintf('Preparing %d remaining scenarios in Batch %d for parallel execution...\n', numMissing, targetBatch);
            
            simIn = repmat(Simulink.SimulationInput('hopper_6dof_NED_v2'), numMissing, 1);
            for k = 1:numMissing
                localIdx = missingLocalIdx(k);
                globalIdx = globalIndices(localIdx);
                currentScenario = scenarioStructs(localIdx);
                
                simIn(k) = simIn(k).setVariable('currentScenario', currentScenario);
                simIn(k) = simIn(k).setPreSimFcn(@(in) local_pre_sim(in, currentScenario));
                simIn(k) = simIn(k).setPostSimFcn(@(out) local_post_sim(out, currentScenario, globalIdx, ckptDir));
            end
            
            fprintf('Running parsim across workers for Batch %d...\n', targetBatch);
            simOuts = parsim(simIn, 'ShowProgress', 'on');
            
            % Reload updated individual checkpoint files into resultsCell
            for k = 1:numMissing
                localIdx = missingLocalIdx(k);
                globalIdx = globalIndices(localIdx);
                ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', globalIdx));
                if exist(ckptFile, 'file')
                    data = load(ckptFile, 'localResult');
                    resultsCell{localIdx} = data.localResult;
                end
            end
        end
        
        elapsedTime = toc;
        fprintf('Completed parallel execution phase for Batch %d in %.2f seconds.\n', targetBatch, elapsedTime);
        
        % --- Export Batch Results ---
        batchResults = [resultsCell{:}];
        batchResultsFile = fullfile(batchDir, sprintf('batch_%03d_results.mat', targetBatch));
        save(batchResultsFile, 'batchResults', 'globalIndices');
        fprintf('Batch results successfully saved to %s\n', batchResultsFile);
        
    case 'serial'
        tic;
        for k = 1:batchNumScen
            if ~isempty(resultsCell{k})
                continue;
            end
            
            globalIdx = globalIndices(k);
            currentScenario = scenarioStructs(k);
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
            
            resultsCell{k} = localResult;
            ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', globalIdx));
            save(ckptFile, 'localResult');
        end
        
        elapsedTime = toc;
        fprintf('Completed sequential execution for Batch %d in %.2f seconds.\n', targetBatch, elapsedTime);
        
        batchResults = [resultsCell{:}];
        batchResultsFile = fullfile(batchDir, sprintf('batch_%03d_results.mat', targetBatch));
        save(batchResultsFile, 'batchResults', 'globalIndices');
        fprintf('Batch results successfully saved to %s\n', batchResultsFile);
        
   case 'nominal'
        tic;
        try
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

function out = local_post_sim(out, scenario, globalIdx, ckptDir)
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
    
    ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', globalIdx));
    save(ckptFile, 'localResult');
end
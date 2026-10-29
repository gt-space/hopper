% --- CTV Monte Carlo Simulation Master Script (Per-Run Checkpointed) ---
% Port of sims/hopper_6dof_sizer/mc_parallelization.m for the CTV sim.
%
% n and runMode are picked up from the caller's workspace if they are already
% set, so this runs unattended under matlab -batch:
%   runMode = 'serial'; n = 20; mc_parallelization_ctv
% Options for runMode: 'serial_isolated', 'parallel', 'serial', 'nominal'.
%
% Use 'serial_isolated' for real Monte Carlo runs. hopper_6dof_NED_v2_ctv is
% only repeatable on the FIRST sim of a MATLAB process: later sims in the same
% process drift by several m/s in landing velocity for bit-identical inputs,
% and no reset fixes it (base workspace wipe, bdclose, clear
% functions/persistent/mex, FastRestart off and Simscape artifact caching off
% were all tried). 'serial_isolated' runs every scenario as sim #1 of a fresh
% process, which costs MATLAB startup per run but is reproducible. 'serial'
% is kept for quick smoke tests only; its results carry an ordinal-dependent
% offset.
warning('off', 'MATLAB:Python:PyNotFound')

% --- Setup CTV-only paths ---
% NOT addpath(genpath(repoRoot)): sims/hopper_6dof_sizer holds same-named
% files (mission_inputs, Outputs, lateralDampingCalc, cg_moi_init, ...) that
% would shadow CTV's.
currentDir = fileparts(mfilename('fullpath'));
ctvPaths = ctv_path_list(currentDir);
addpath(ctvPaths{:});
fprintf('mc_parallelization_ctv path check:\n');
fprintf('  mission_inputs     -> %s\n', which('mission_inputs'));
fprintf('  lateralDampingCalc -> %s\n', which('lateralDampingCalc'));

mdl = 'hopper_6dof_NED_v2_ctv';

if ~exist('runMode', 'var') || isempty(runMode)
    runMode = 'parallel';
end
if ~exist('n', 'var') || isempty(n)
    n = 100;
end
if ~strcmp(runMode, 'nominal') && (~isnumeric(n) || ~isscalar(n) || n <= 0)
    error('n must be a positive integer number of Monte Carlo scenarios.');
end
fprintf('Model: %s | runMode: %s | n: %d\n', mdl, runMode, n);

% The real work runs in a local function, not here. mc_sim_setup_ctv has to
% run LinerizationMaster_ctv and Trajectory2 in the base workspace (they are
% scripts that read and write it), and those two assign i, t, m, x, z, A, B,
% scenario and more. A driver loop living in the base workspace would have its
% own loop index overwritten mid-iteration; a function workspace is immune.
run_mc(mdl, runMode, n, currentDir, ctvPaths);

function run_mc(mdl, runMode, n, currentDir, ctvPaths)

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

% Attach all required models, lookup files, and subfolder directories to
% workers. CTV keeps its inputs in subfolders, so this search is recursive.
addAllFiles = {};
for pat = {'*.m', '*.mat', '*.slx', '*.xlsx', '*.json'}
    f = dir(fullfile(currentDir, '**', pat{1}));
    addAllFiles = [addAllFiles, ...
        cellfun(@(d, nm) fullfile(d, nm), {f.folder}, {f.name}, ...
                'UniformOutput', false)];  %#ok<AGROW>
end
if strcmp(runMode, 'parallel')
    p = gcp();
    addAttachedFiles(p, addAllFiles);
    % The CTV subfolders have to be on each worker's path before the model
    % loads, which happens before the PreSimFcn runs.
    parfevalOnAll(@addpath, 0, ctvPaths{:});
end

jsonFile = fullfile(currentDir, 'mc_params_ctv.json');
ckptDir = fullfile(currentDir, 'mc_checkpoints'); % Use absolute path

if ~exist(ckptDir, 'dir')
    mkdir(ckptDir);
end

if ~strcmp(runMode, 'nominal')
    [scenarios, mcTable] = generateScenarios_ctv(jsonFile, n, currentDir);
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

            simIn = repmat(Simulink.SimulationInput(mdl), numMissing, 1);
            for k = 1:numMissing
                i = missingIdx(k);
                currentScenario = scenarioStructs(i);

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
        nominalScenario = nominal_scenario(jsonFile);
        mc_sim_setup_ctv(nominalScenario);
        simInputNom = Simulink.SimulationInput(mdl);
        simOutNom = sim(simInputNom);
        nominal = apply_constraints(mc_main_ctv(nominalScenario, simOutNom));

        % --- Export Final Results ---
        resultsCell = cellfun(@apply_constraints, resultsCell, 'UniformOutput', false);
        results = [resultsCell{:}];
        save(fullfile(currentDir, 'mc_results_ctv_parallel.mat'), 'results', 'nominal', 'scenarios');
        summarize(results);
        fprintf('Final results and nominal data saved to mc_results_ctv_parallel.mat\n');

    case 'serial_isolated'
        tic;
        scenFile = fullfile(ckptDir, 'scenarios.mat');
        nominalScenario = nominal_scenario(jsonFile);
        save(scenFile, 'scenarioStructs', 'nominalScenario');
        matlabExe = fullfile(matlabroot, 'bin', 'matlab');

        missingIdx = find(cellfun(@isempty, resultsCell));
        fprintf('Running %d scenarios, one MATLAB process each...\n', numel(missingIdx));
        for k = 1:numel(missingIdx)
            idx = missingIdx(k);
            tRun = tic;
            run_isolated(idx, scenFile, ckptDir, currentDir, matlabExe);
            ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', idx));
            if exist(ckptFile, 'file')
                data = load(ckptFile, 'localResult');
                resultsCell{idx} = data.localResult;
                tag = 'ok';
                if ~data.localResult.status.success
                    tag = ['ERROR: ' data.localResult.status.error];
                end
            else
                tag = 'NO CHECKPOINT (worker died)';
            end
            fprintf('  %d/%d index %d: %s in %.1f s\n', ...
                k, numel(missingIdx), idx, tag, toc(tRun));
        end

        elapsedTime = toc;
        fprintf('Completed isolated execution in %.2f seconds.\n', elapsedTime);

        fprintf('Running nominal case...\n');
        run_isolated(0, scenFile, ckptDir, currentDir, matlabExe);
        nomFile = fullfile(ckptDir, 'nominal.mat');
        if ~exist(nomFile, 'file')
            error('Nominal run produced no checkpoint.');
        end
        nominal = apply_constraints(load(nomFile, 'localResult').localResult);

        resultsCell = cellfun(@apply_constraints, resultsCell, 'UniformOutput', false);
        results = [resultsCell{:}];
        save(fullfile(currentDir, 'mc_results_ctv_serial_isolated.mat'), ...
            'results', 'nominal', 'scenarios');
        summarize(results);
        fprintf('Final results saved to mc_results_ctv_serial_isolated.mat\n');

    case 'serial'
        tic;
        for i = 1:n
            if ~isempty(resultsCell{i})
                continue;
            end

            currentScenario = scenarioStructs(i);
            localResult = struct();

            try
                mc_sim_setup_ctv(currentScenario);
                simInput = Simulink.SimulationInput(mdl);

                sim_out = sim(simInput);
                localResult = mc_main_ctv(currentScenario, sim_out);
            catch ME
                localResult.scenario = currentScenario;
                localResult.status.success = false;
                localResult.status.error = ME.message;
                localResult.status.pass = false;
            end
            localResult = apply_constraints(localResult);

            resultsCell{i} = localResult;

            % Save individual checkpoint file immediately after each run
            ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', i));
            save(ckptFile, 'localResult');
        end

        elapsedTime = toc;
        fprintf('Completed sequential execution runs in %.2f seconds.\n', elapsedTime);

        % --- Run Nominal Case for Comparison ---
        fprintf('Running nominal case...\n');
        nominalScenario = nominal_scenario(jsonFile);
        mc_sim_setup_ctv(nominalScenario);
        simOutNom = sim(Simulink.SimulationInput(mdl));
        nominal = apply_constraints(mc_main_ctv(nominalScenario, simOutNom));

        resultsCell = cellfun(@apply_constraints, resultsCell, 'UniformOutput', false);
        results = [resultsCell{:}];
        save(fullfile(currentDir, 'mc_results_ctv_serial.mat'), 'results', 'nominal', 'scenarios');
        summarize(results);
        fprintf('Final results successfully saved to mc_results_ctv_serial.mat\n');

   case 'nominal'
        tic;
        nominalScenario = nominal_scenario(jsonFile);

        mc_sim_setup_ctv(nominalScenario);

        simInput = Simulink.SimulationInput(mdl);
        sim_out = sim(simInput);
        nominal = apply_constraints(mc_main_ctv(nominalScenario, sim_out));

        save(fullfile(currentDir, 'mc_results_ctv_nominal.mat'), 'nominal');
        disp('Nominal simulation completed and saved successfully.');

        elapsedTime = toc;
        fprintf('Nominal run completed in %.2f seconds.\n', elapsedTime);

    otherwise
        error(['Invalid runMode specified. Use ''serial_isolated'', ' ...
            '''parallel'', ''serial'' or ''nominal''.']);
end

end

function run_isolated(idx, scenFile, ckptDir, ctvDir, matlabExe)
    % One scenario in its own MATLAB process. The model prints tens of
    % thousands of lines per run, so the worker's output goes to a log that is
    % deleted unless the run fails.
    logFile = fullfile(ckptDir, sprintf('worker_%05d.log', idx));
    cmd = sprintf('"%s" -batch "addpath(''%s''); mc_run_one_ctv(%d, ''%s'', ''%s'')" > "%s" 2>&1', ...
        matlabExe, ctvDir, idx, scenFile, ckptDir, logFile);
    st = system(cmd);
    if st == 0 && exist(logFile, 'file')
        delete(logFile);
    elseif st ~= 0
        fprintf('    worker exited %d; log kept at %s\n', st, logFile);
    end
end

function paths = ctv_path_list(ctvDir)
    paths = [{ctvDir}, cellfun(@(s) fullfile(ctvDir, s), ...
        {'sizing', 'inputs', 'propulsion', 'dynamics', 'loads', 'cg', 'gnc', 'baffle'}, ...
        'UniformOutput', false)];
end

function scenario = nominal_scenario(jsonFile)
    % Nominal values straight out of the bounds file, so the nominal run uses
    % the same parameter set as the dispersed ones.
    bounds = readstruct(jsonFile);
    names = fieldnames(bounds);
    scenario = struct();
    for k = 1:length(names)
        scenario.(names{k}) = bounds.(names{k}).nominal;
    end
end

function in = local_pre_sim(in, scenario)
    mc_sim_setup_ctv(scenario);
end

function out = local_post_sim(out, scenario, idx, ckptDir)
    % This runs on the worker thread the exact moment a simulation finishes
    try
        if ~isempty(out.ErrorMessage)
            error(out.ErrorMessage);
        end
        localResult = mc_main_ctv(scenario, out);
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

function r = apply_constraints(r)
    % Guarantee every result carries status.pass, including the error cases
    % that never got far enough to produce vehicle/flight fields.
    if isempty(r)
        return
    end
    if isfield(r, 'status') && isfield(r.status, 'success') && r.status.success
        r = check_constraints(r);
    else
        if ~isfield(r, 'status') || ~isfield(r.status, 'pass')
            r.status.pass = false;
        end
        r.status.pass_flight = false;
        r.status.failed = {'errored'};
    end
end

function result = check_constraints(result)
    % Criteria from sims/hopper_6dof_sizer/mc_runner.m, unchanged. status.pass
    % keeps the original meaning; status.failed names the criteria that broke
    % and status.pass_flight drops the two TWR criteria, which CTV cannot meet
    % at a final TWR of 0.74 and so would otherwise mask every flight result.
    v = result.vehicle;
    f = result.flight;

    % name, violated, counts_towards_pass_flight
    checks = { ...
        'twr_initial',  v.twr_initial     < 1.0, false; ...
        'twr_final',    v.twr_final       < 1.0, false; ...
        'dry_mass',     v.dry_mass        <= 0,  true;  ...
        'wet_mass',     v.wet_mass        <= 0,  true;  ...
        'max_altitude', f.max_altitude    < 51,  true;  ...
        'max_ascent_vel',  f.max_ascent_vel  > 7,   true;  ...
        'max_descent_vel', f.max_descent_vel > 7,   true;  ...
        'landing_vel',  f.landing_vel     > 0.5, true};

    violated = cell2mat(checks(:, 2));
    inFlight = cell2mat(checks(:, 3));

    result.status.pass        = ~any(violated);
    result.status.pass_flight = ~any(violated & inFlight);
    result.status.failed      = reshape(checks(violated, 1), 1, []);
end

function summarize(results)
    if isempty(results); return; end
    nr = numel(results);
    ok = arrayfun(@(r) isfield(r, 'status') && r.status.success, results);
    pass = arrayfun(@(r) isfield(r, 'status') && isfield(r.status, 'pass') && r.status.pass, results);
    passF = arrayfun(@(r) isfield(r, 'status') && isfield(r.status, 'pass_flight') && r.status.pass_flight, results);
    fprintf('\n%s\n', repmat('=', 1, 52));
    fprintf('Scenarios:    %d\n', nr);
    fprintf('Pass:         %d (%.1f%%)\n', sum(pass), 100*sum(pass)/nr);
    fprintf('Fail:         %d (%.1f%%)\n', sum(ok & ~pass), 100*sum(ok & ~pass)/nr);
    fprintf('Pass(flight): %d (%.1f%%)   [TWR criteria excluded]\n', sum(passF), 100*sum(passF)/nr);
    fprintf('Errors:       %d (%.1f%%)\n', sum(~ok), 100*sum(~ok)/nr);

    % per-constraint tally
    allNames = {};
    for i = 1:nr
        if isfield(results(i), 'status') && isfield(results(i).status, 'failed')
            allNames = [allNames, results(i).status.failed];  %#ok<AGROW>
        end
    end
    if ~isempty(allNames)
        u = unique(allNames);
        fprintf('%s\n', repmat('-', 1, 52));
        for i = 1:numel(u)
            c = sum(strcmp(allNames, u{i}));
            fprintf('  broke %-18s %4d (%.1f%%)\n', u{i}, c, 100*c/nr);
        end
    end
    fprintf('%s\n', repmat('=', 1, 52));
end

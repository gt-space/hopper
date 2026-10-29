function mc_run_one_ctv(idx, scenFile, ckptDir)
%MC_RUN_ONE_CTV Run one CTV Monte Carlo scenario and checkpoint the result.
%   Entry point for the 'serial_isolated' mode of mc_parallelization_ctv,
%   which launches one MATLAB process per scenario.
%
%   idx      scenario index into scenarioStructs in scenFile, or 0 for the
%            nominal scenario
%   scenFile .mat holding scenarioStructs and nominalScenario
%   ckptDir  directory to write scenario_%05d.mat (or nominal.mat) into
%
%   hopper_6dof_NED_v2_ctv is only repeatable on the first sim of a MATLAB
%   process: later sims in the same process drift by several m/s in landing
%   velocity for bit-identical inputs. Running every scenario as sim #1 of a
%   fresh process is what makes the Monte Carlo reproducible.

ctvDir = fileparts(mfilename('fullpath'));
addpath(ctvDir);
mdl = 'hopper_6dof_NED_v2_ctv';

S = load(scenFile);
if idx == 0
    scenario = S.nominalScenario;
    ckptFile = fullfile(ckptDir, 'nominal.mat');
else
    scenario = S.scenarioStructs(idx);
    ckptFile = fullfile(ckptDir, sprintf('scenario_%05d.mat', idx));
end

localResult = struct();
try
    mc_sim_setup_ctv(scenario);
    sim_out = sim(Simulink.SimulationInput(mdl));
    localResult = mc_main_ctv(scenario, sim_out);
catch ME
    localResult.scenario = scenario;
    localResult.status.success = false;
    localResult.status.error = ME.message;
    localResult.status.pass = false;
end

save(ckptFile, 'localResult');
fprintf('mc_run_one_ctv: wrote %s (success=%d)\n', ckptFile, ...
    localResult.status.success);
end

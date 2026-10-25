function W = hopperWorkspace(scenario)
%HOPPERWORKSPACE  Every variable the hopper setup creates, as a struct.
%   W = hopperWorkspace()          nominal vehicle (what sim_setup builds)
%   W = hopperWorkspace(scenario)  one Monte Carlo scenario (mc_sim_setup)
%
%   The setup scripts (main, load_lookup, cg_moi_test, LinerizationMaster,
%   Trajectory2, ...) still talk to each other through the base workspace.
%   This function runs them there in isolation: the caller's base variables
%   are saved first and put back afterwards, so calling it has no side
%   effects. Use buildSimInput to turn W into a Simulink.SimulationInput.

saved = snapshotBase();
restore = onCleanup(@() restoreBase(saved));
evalin('base', 'clearvars');

if nargin < 1 || isempty(scenario)
    evalin('base', 'sim_setup');
else
    mc_sim_setup(scenario);
end
W = snapshotBase();
end

function S = snapshotBase()
S = struct();
for n = evalin('base', 'who')'
    S.(n{1}) = evalin('base', n{1});
end
end

function restoreBase(S)
evalin('base', 'clearvars');
for n = fieldnames(S)'
    assignin('base', n{1}, S.(n{1}));
end
end

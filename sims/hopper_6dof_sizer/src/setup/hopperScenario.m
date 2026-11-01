function s = hopperScenario(set)
%HOPPERSCENARIO  The Monte Carlo scenario main.m passes to mission_inputs.
%   s = hopperScenario()          current scenario ([] = nominal)
%   old = hopperScenario(scen)    set one, returning the previous one
%
%   Lets sim_setup (which starts with clearvars) build one dispersed run.
%   hopperWorkspace uses it for vehicles whose Monte Carlo goes through
%   sim_setup (the CTV); the hopper's still goes through mc_sim_setup.

persistent current
s = current;
if nargin > 0
    current = set;
end
end

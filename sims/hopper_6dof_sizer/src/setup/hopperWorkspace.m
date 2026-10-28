function W = hopperWorkspace(scenario, variants)
%HOPPERWORKSPACE  Every variable the hopper setup creates, as a struct.
%   W = hopperWorkspace()          nominal vehicle (what sim_setup builds)
%   W = hopperWorkspace(scenario)  one Monte Carlo scenario (mc_sim_setup)
%   W = hopperWorkspace(scenario, variants)
%                                  for a variantConfig: VEHICLE picks the
%                                  parameter set (hopperVehicle), and
%                                  CONTROLLER=lqi adds the LQI gains
%
%   The setup scripts (main, load_lookup, cg_moi_test, LinerizationMaster,
%   Trajectory2, ...) still talk to each other through the base workspace.
%   This function runs them there in isolation: the caller's base variables
%   are saved first and put back afterwards, so calling it has no side
%   effects. Use buildSimInput to turn W into a Simulink.SimulationInput.

if nargin < 2 || isempty(variants), variants = variantConfig(); end
A = variantAxes();
vehicle = A.VEHICLE(variants.VARIANT_VEHICLE);
if vehicle ~= "hopper" && nargin >= 1 && ~isempty(scenario)
    error('hopperWorkspace:mcVehicle', ...
        'Monte Carlo scenarios (mc_sim_setup) only model the hopper so far, not the %s.', vehicle);
end

saved = snapshotBase();
restore = onCleanup(@() restoreBase(saved));
oldVehicle = hopperVehicle(vehicle);
restoreVehicle = onCleanup(@() hopperVehicle(oldVehicle));
evalin('base', 'clearvars');

if nargin < 1 || isempty(scenario)
    evalin('base', 'sim_setup');
else
    mc_sim_setup(scenario);
end
W = snapshotBase();

if A.CONTROLLER(variants.VARIANT_CONTROLLER) == "lqi"
    [W.Klqi, W.p] = lqiGains();
end

% Data a model callback loads into base itself: the fluid network's InitFcn
% loads FuelProps (Simscape tank fluid properties). buildSimInput only sets
% what a model reads, so models without the fluid network are unaffected.
F = load('FuelPropsData.mat');
for f = fieldnames(F)'
    W.(f{1}) = F.(f{1});
end
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

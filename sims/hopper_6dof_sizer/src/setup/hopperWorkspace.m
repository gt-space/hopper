function W = hopperWorkspace(scenario, variants)
%HOPPERWORKSPACE  Every variable the hopper setup creates, as a struct.
%   W = hopperWorkspace()          nominal vehicle (what sim_setup builds)
%   W = hopperWorkspace(scenario)  one Monte Carlo scenario (mc_sim_setup)
%   W = hopperWorkspace(scenario, variants)
%                                  for a variantConfig: VEHICLE picks the
%                                  parameter set (hopperVehicle), and
%                                  CONTROLLER=lqi adds the LQI gains. A CTV
%                                  scenario runs sim_setup with its
%                                  dispersions and the nominal CTV's
%                                  reference and gains
%
%   The setup scripts (main, load_lookup, cg_moi_test, LinerizationMaster,
%   Trajectory2, ...) still talk to each other through the base workspace.
%   This function runs them there in isolation: the caller's base variables
%   are saved first and put back afterwards, so calling it has no side
%   effects. Use buildSimInput to turn W into a Simulink.SimulationInput.

if nargin < 1, scenario = []; end
if nargin < 2 || isempty(variants), variants = variantConfig(); end
A = variantAxes();
vehicle = A.VEHICLE(variants.VARIANT_VEHICLE);

saved = snapshotBase();
restore = onCleanup(@() restoreBase(saved));
oldVehicle = hopperVehicle(vehicle);
restoreVehicle = onCleanup(@() hopperVehicle(oldVehicle));
evalin('base', 'clearvars');

if isempty(scenario)
    evalin('base', 'sim_setup');
elseif vehicle == "hopper"
    mc_sim_setup(scenario);
else
    % Other vehicles: sim_setup with the scenario's dispersions (their
    % inputs file applies them), flying the controller designed for the
    % nominal vehicle, as mc_sim_setup does for the hopper
    W = withScenario(scenario);
    C = nominalController(vehicle);
    for f = fieldnames(C)'
        W.(f{1}) = C.(f{1});
    end
end
if isempty(scenario) || vehicle == "hopper"
    W = snapshotBase();
end

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

function W = withScenario(scenario)
old = hopperScenario(scenario);
restoreScenario = onCleanup(@() hopperScenario(old));
evalin('base', 'sim_setup');
W = snapshotBase();
end

function C = nominalController(vehicle)
% Reference trajectory and gains for the nominal vehicle, cached per
% session (and per parsim worker) since every scenario shares them
persistent cache
if isempty(cache), cache = containers.Map(); end
key = char(vehicle);
if ~isKey(cache, key)
    evalin('base', 'clearvars');
    evalin('base', 'sim_setup');
    W = snapshotBase();
    names = {'tgrid', 'unom', 'Target_Trajectory', 'K1flat', 'K2grid', 'K2grid_mod', ...
             'z_ref_ts', 'vz_ref_ts', 'T_ts', 'T_ref_ts'};
    C = struct();
    for n = names(isfield(W, names))
        C.(n{1}) = W.(n{1});
    end
    cache(key) = C;
end
C = cache(key);
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

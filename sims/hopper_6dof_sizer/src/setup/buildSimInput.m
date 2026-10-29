function [in, W] = buildSimInput(model, scenario, W, variants)
%BUILDSIMINPUT  Simulink.SimulationInput with every variable a model reads.
%   in = buildSimInput()                    nominal run of hopper_6dof_NED_v2
%   in = buildSimInput(model)               nominal run of another model
%   in = buildSimInput(model, scenario)     one Monte Carlo scenario
%   in = buildSimInput(model, scenario, W)  reuse a hopperWorkspace result
%                                           (add extra variables to W first)
%   in = buildSimInput(model, scenario, W, variants)
%                                           variant choices for
%                                           models/hopper_6dof.slx, from
%                                           variantConfig (default: all
%                                           defaults, same as v2)
%   [in, W] = buildSimInput(...)            also return the workspace struct
%
%   With FDIR on, StopTime is the trajectory end + 30 s so an abort to
%   land can finish.
%
%   model can also be an existing Simulink.SimulationInput (for example
%   inside a parsim PreSimFcn); the variables are added to it.
%
%   The variables are set on the SimulationInput, so out = sim(in) needs
%   nothing in the base workspace and parsim can run inputs side by side.
%   Which variables a model reads is found once per session with
%   Simulink.findVars and cached (see modelVariables).

if nargin < 1 || isempty(model), model = 'hopper_6dof_NED_v2'; end
if nargin < 2, scenario = []; end
if nargin < 4 || isempty(variants), variants = variantConfig(); end
if nargin < 3 || isempty(W), W = hopperWorkspace(scenario, variants); end
for f = fieldnames(variants)'
    W.(f{1}) = variants.(f{1});   % only models with variants read these
end
if variants.VARIANT_FDIR == 2 || variants.VARIANT_NAVIGATION > 1 || variants.VARIANT_FAULT_INJECTION == 2
    W = sensorWorkspace(W, variants); % the shared Sensors block: FDIR and navigation read it
end                                   % (keeps a fault table already in W.faults)
if variants.VARIANT_NAVIGATION > 1
    W = navWorkspace(W, variants);
end
if variants.VARIANT_FDIR == 2 || variants.VARIANT_FAULT_INJECTION == 2
    W = fdirWorkspace(W, variants);   % reads the sensor and nav parameters above
end

if isa(model, 'Simulink.SimulationInput')
    in = model;
    model = in.ModelName;
else
    in = Simulink.SimulationInput(model);
end
for v = modelVariables(model, W)
    if ~isfield(W, v{1})
        error('buildSimInput:missingVariable', ...
            '%s reads ''%s'', which the setup did not create.', model, v{1});
    end
    in = in.setVariable(v{1}, W.(v{1}));
end
if variants.VARIANT_FDIR == 2 && isfield(W, 'tgrid')
    % An abort to land can start late and high; leave it time to come
    % down (the touchdown stop still ends a normal flight on time)
    in = in.setModelParameter('StopTime', num2str(W.tgrid(end) + 30));
end
end

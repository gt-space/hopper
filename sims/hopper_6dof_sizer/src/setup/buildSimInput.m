function [in, W] = buildSimInput(model, scenario, W)
%BUILDSIMINPUT  Simulink.SimulationInput with every variable a model reads.
%   in = buildSimInput()                    nominal run of hopper_6dof_NED_v2
%   in = buildSimInput(model)               nominal run of another model
%   in = buildSimInput(model, scenario)     one Monte Carlo scenario
%   in = buildSimInput(model, scenario, W)  reuse a hopperWorkspace result
%                                           (add extra variables to W first)
%   [in, W] = buildSimInput(...)            also return the workspace struct
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
if nargin < 3 || isempty(W), W = hopperWorkspace(scenario); end

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
end

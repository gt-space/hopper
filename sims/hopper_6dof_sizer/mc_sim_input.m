function in = mc_sim_input(model, flags)
%MC_SIM_INPUT  SimulationInput for the chosen 6DOF model, used by the Monte Carlo scripts.
%   in = mc_sim_input(model)          full physics
%   in = mc_sim_input(model, flags)   flags = [enable_slosh enable_wind], toggleable only
%
%   model:
%     'original'   hopper_6dof_NED_v2             full physics (the master)
%     'toggleable' hopper_6dof_NED_v2_simplified  slosh / wind switchable with flags
%     'minimal'    hopper_6dof_NED_v2_minimal     no slosh or wind, constant mass;
%                                                 does NOT apply cg_factor dispersions
%
%   The toggles are set on the SimulationInput, so they apply to that one run and
%   never change the base workspace.

if nargin < 1 || isempty(model), model = 'original'; end
if nargin < 2 || isempty(flags), flags = [1 1];      end

switch lower(model)
    case 'original',   name = 'hopper_6dof_NED_v2';
    case 'toggleable', name = 'hopper_6dof_NED_v2_simplified';
    case 'minimal',    name = 'hopper_6dof_NED_v2_minimal';
    otherwise
        error('mc_sim_input:unknownModel', ...
            'Unknown model "%s". Use ''original'', ''toggleable'' or ''minimal''.', model);
end

in = Simulink.SimulationInput(name);
if strcmpi(model, 'toggleable')
    in = in.setVariable('enable_slosh', flags(1));
    in = in.setVariable('enable_wind',  flags(2));
elseif ~isequal(flags(:)', [1 1])
    warning('mc_sim_input:flagsIgnored', ...
        'Flags only apply to the toggleable model; ignoring them for "%s".', model);
end
end

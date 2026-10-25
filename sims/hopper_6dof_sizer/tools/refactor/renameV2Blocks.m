function renamed = renameV2Blocks(model)
%RENAMEV2BLOCKS  Give v2's generically named blocks names that say what they do.
%   renamed = renameV2Blocks()        renames in hopper_6dof_NED_v2 and saves
%   renamed = renameV2Blocks(model)   any model built from v2 (for example
%                                     the fluid-integration model); paths
%                                     the model doesn't have are skipped
%
%   Only names change; wiring and parameters don't. Each name was read off
%   the block's code or logic (stage 5, step 1). Scripts that address blocks
%   by path (validation/mass_properties/massPropertiesCheck.m,
%   tools/refactor/build_v2_fluid_integration.m) use the new names.
%   Returns {old path, new name} for every block renamed.

if nargin < 1, model = 'hopper_6dof_NED_v2'; end

% {path relative to the model (with '//' for a '/' in a name), new name}
top = {
    'Subsystem',                       'Thrust_to_Mdot_LUT'           % thrust -> ox/fuel mdot lookups
    'Subsystem1',                      'Propellant_Mass_Integration'  % integrate mdot -> ox/fuel mass
    'Subsystem2',                      'Propellant_Depletion_Cutoff'  % thrust = 0 once either tank is empty
    'Subsystem3',                      'Touchdown_Stop'               % stop the sim at ground contact, descending
    'Subsystem4',                      'Landing_Engine_Cutoff'        % engine cutoff just above the ground
    'MATLAB Function',                 'LQR_Gain_Scheduler'           % Gain_Selector: time-scheduled LQR
    'CG // MOI Script v3',             'Mass_Properties'              % CG, inertia, dI/dt (+ slosh diagnostics)
    sprintf('Control Inputs -> Fxyz & Mxyz'), 'Actuators_TVC_RCS'     % TVC/RCS actuator models -> body F, M
    'Fluid Network2',                  'Fluid_Network'                % Simscape feed system + engine (off in v2)
    'Ground ',                         'Ground_Contact_Unwired'       % old leg model; outputs only reach scopes
    };
inner = {
    'Mass_Properties/MATLAB Function',   'Ox_Slosh_Inertia'           % zeroed: propellant inertia is in the LUTs
    'Mass_Properties/MATLAB Function1',  'Fu_Slosh_Inertia'
    'Mass_Properties/MATLAB Function2',  'TVC_Skew_Matrix'            % [0 -dz dy; dz 0 0; -dy 0 0]
    'Fluid_Network/MATLAB Function',     'Engine_Model'               % mdot, c*, Cf -> thrust, Pc, MR
    'Fluid_Network/Fuel Throttle Valve Controller/MATLAB Function5', 'Valve_Rate_Limit_Schedule'
    'Fluid_Network/Ox Throttle Valve Controller/MATLAB Function5',   'Valve_Rate_Limit_Schedule'
    'Fluid_Network/Thrust Commanded to mdot + bunch of other bs',    'Thrust_Cmd_to_Mdot'
    'Slosh/MATLAB Function',             'Ox_Slosh_Parameters'        % SP-106 mode frequency, mass, length
    'Slosh/MATLAB Function1',            'Fu_Slosh_Parameters'
    'Slosh/Subsystem',                   'Ox_Effective_Accel'         % a_cg + rotational terms at the tank
    'Slosh/Subsystem1',                  'Fu_Effective_Accel'
    sprintf('Propellant_Depletion_Cutoff/If Action\nSubsystem'),  'Has_Propellant'
    sprintf('Propellant_Depletion_Cutoff/If Action\nSubsystem1'), 'Out_Of_Propellant'
    };
% The six pendulum blocks share one layout
for tank = {'Fuel', 'Oxidizer'}
    for axis = {'Axial', 'Pitch', 'Yaw'}
        base = sprintf('Slosh/%s %s Force/', tank{1}, axis{1});
        inner(end+1:end+5, :) = {
            [base 'MATLAB Function'],  'Accel_Term'                   % -a / l
            [base 'MATLAB Function1'], 'Stiffness_Term'               % -w^2 phi
            [base 'MATLAB Function2'], 'Damping_Term'                 % -2 zeta w phi_dot
            [base 'MATLAB Function3'], 'Pendulum_Force'               % -m l phi_ddot on the tank
            [base 'MATLAB Function4'], 'Pendulum_Torque'              % zeroed: hinge carries no moment
            };
    end
end

load_system(model);
renamed = {};
for M = {top, inner}
    map = M{1};
    for k = 1:size(map, 1)
        old = [model '/' map{k, 1}];
        if getSimulinkBlockHandle(old) < 0, continue; end
        set_param(old, 'Name', map{k, 2});
        renamed(end+1, :) = {map{k, 1}, map{k, 2}}; %#ok<AGROW>
    end
end
if ~isempty(renamed)
    save_system(model);
end
fprintf('%s: %d blocks renamed\n', model, size(renamed, 1));
end

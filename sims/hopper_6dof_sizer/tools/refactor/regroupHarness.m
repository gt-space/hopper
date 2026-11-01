function regroupHarness(model)
%REGROUPHARNESS  The harness top level as seven groups joined by named buses.
%   regroupHarness()        models/hopper_6dof.slx (saves it)
%   regroupHarness(model)   another loaded copy (saves it)
%
%   The top level had grown to ~96 blocks plus 43 commented-out leftovers.
%   This deletes the leftovers and moves every other block, unchanged, into
%   a virtual subsystem by area (moveIntoGroup):
%     Environment         wind forces and moments, random roll disturbance
%     Navigation_Sensors  shared sensor model and the NAVIGATION variant
%     GNC                 controller, FDIR, actuator fault injection
%     Actuation           thrust command chain and TVC/RCS actuators
%     Mass_Propulsion     engine cutoffs, propulsion, propellant, mass properties
%     Logging             every To Workspace, scope and their references
%     Vehicle_Dynamics    6DOF, slosh, ground contact, force/moment sums
%   then joins the groups with one bus each, every signal named
%   (busifyHarness), and lays the result out (layoutHarnessGroups). The
%   default flight is bit-identical before and after. Rerunning is harmless.

if nargin < 1, model = 'hopper_6dof'; end
load_system(model);
if getSimulinkBlockHandle([model '/Vehicle_Dynamics']) < 0
    deleteLeftovers(model);
    G = groups();
    for g = 1:size(G, 1), moveIntoGroup(model, G{g, 1}, G{g, 2}); end
    save_system(model);
else
    fprintf('%s: already grouped\n', model);
end
busifyHarness(model);
layoutHarnessGroups(model);
end

function G = groups()
G = {
    'Environment', {'Wind_Forces', 'MatrixMultiply1', 'MatrixMultiply2', 'Bernoulli Binary Generator', 'Random Source', 'Product1'}
    'Navigation_Sensors', {'Sensors', 'Navigation', 'Mux'}
    'GNC', {'Controller', 'FDIR', 'Actuator_Faults'}
    'Actuation', {'Demux', 'Actuators_TVC_RCS', 'IC', 'Rate Limiter1', 'Saturation2', 'Transport Delay2', 'Quantizer', 'Product'}
    'Mass_Propulsion', {'Landing_Engine_Cutoff', 'Propellant_Depletion_Cutoff', 'Propulsion', 'Propellant_Flow', ...
        'Propellant_Mass_Integration', 'Mass_Properties', 'Sum', 'Integrator', 'Sum1', 'Mass', 'Saturation1'}
    'Logging', {'To Workspace', 'To Workspace1', 'To Workspace2', 'To Workspace3', 'To Workspace4', 'To Workspace5', ...
        'To Workspace6', 'To Workspace7', 'To Workspace8', 'To Workspace9', 'To Workspace10', 'To Workspace11', ...
        'To Workspace MOI_with_slosh', 'To Workspace cg_with_slosh', 'CG Accel', 'Pitch', 'Scope6', 'Vz', 'X', 'Y', 'Z', ...
        'From Workspace', 'From Workspace1', 'From Workspace2', 'Radians to Degrees', 'Radians to Degrees1', 'Radians to Degrees2'}
    'Vehicle_Dynamics', {'6DOF (Quaternion)', 'Slosh', 'Ground_Contact', 'Slosh Specific Force', 'Sum2', 'Sum3', 'Sum4', ...
        'Sum5', 'Sum6', 'Sum7', 'Sum8', 'Sum9', 'Sum10', 'Sum11', 'Sum12', 'Sum13', 'Sum14', 'Sum15', 'Demux1', 'Demux2', ...
        'Demux3', 'Demux4', 'Demux5', 'Demux6', 'Demux7', 'Demux10', 'Mux1', 'Mux2', 'Cross Product', 'Divide', ...
        'Unit Delay', 'Unit Delay1', 'Derivative3', 'MatrixMultiply', 'g', 'Invert  3x3 Matrix', 'Touchdown_Stop'}
    };
end

function deleteLeftovers(model)
% commented-out blocks, and IC1, whose output fed only commented blocks
top = find_system(model, 'SearchDepth', 1, 'IncludeCommented', 'on', 'Type', 'Block');
top = top(~strcmp(top, model));
dead = top(~strcmp(get_param(top, 'Commented'), 'off'));
if getSimulinkBlockHandle([model '/IC1']) > 0, dead{end+1} = [model '/IC1']; end
for k = 1:numel(dead), delete_block(dead{k}); end
fprintf('%s: %d leftover blocks deleted\n', model, numel(dead));
end

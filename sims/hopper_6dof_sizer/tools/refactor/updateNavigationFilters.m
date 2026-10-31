function updateNavigationFilters()
%UPDATENAVIGATIONFILTERS  Filter-model changes the harness needs (nav team's files).
%   Hopper_GMEKF.slx and CTV_GMEKF.slx (subsystems/navigation/Filters):
%   - the Correct step returns its per-sensor diagnostics (accepted flags
%     and NIS, gmekf_correct_hopper / gmekf_correct) on a third outport,
%     'diagnostics', which FDIR's estimator checks read
%   - CTV_GMEKF uses a fixed-step solver like Hopper_GMEKF: a model
%     referenced from the fixed-step harness must be fixed-step
%   Rerunning is harmless.

filters = {'Hopper_GMEKF', 'gmekf_correct_hopper'; 'CTV_GMEKF', 'gmekf_correct'};
for k = 1:size(filters, 1)
    m = filters{k, 1};
    load_system(m);
    changed = false;
    if ~strcmp(get_param(m, 'SolverType'), 'Fixed-step')
        set_param(m, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepAuto', 'FixedStep', 'auto');
        changed = true;
        fprintf('%s: fixed-step solver\n', m);
    end
    c = [m '/Correct step'];
    ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', c);
    s = ch.Script;
    if ~contains(s, 'diagnostics')
        s = regexprep(s, 'function \[xnp1_c, Pnp1_c\] =', 'function [xnp1_c, Pnp1_c, diagnostics] =', 'once');
        s = regexprep(s, '\[xnp1_c, Pnp1_c\] = ', '[xnp1_c, Pnp1_c, diagnostics] = ', 'once');
        if ~contains(s, ['diagnostics] = ' filters{k, 2}])
            error('updateNavigationFilters:script', '%s: could not find the %s call in the Correct step.', m, filters{k, 2});
        end
        ch.Script = s;
        changed = true;
    end
    if getSimulinkBlockHandle([m '/diagnostics']) < 0
        n = numel(find_system(m, 'SearchDepth', 1, 'BlockType', 'Outport'));
        pos = get_param([m '/covariance'], 'Position');
        add_block('built-in/Outport', [m '/diagnostics'], 'Port', num2str(n + 1), ...
            'Position', pos + [0 60 0 60]);
        add_line(m, 'Correct step/3', 'diagnostics/1', 'autorouting', 'smart');
        changed = true;
        fprintf('%s: diagnostics outport\n', m);
    end
    if changed, save_system(m); end
    close_system(m, 0);
end
end

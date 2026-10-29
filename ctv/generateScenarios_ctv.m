function [scenarios, mcTable] = generateScenarios_ctv(jsonFile, n, outDir)
%GENERATESCENARIOS Latin-hypercube Monte Carlo scenarios from a bounds file.
%   Copy of sims/hopper_6dof_sizer/generateScenarios.m. The only change is
%   outDir: the scenario table is written next to the results rather than to
%   whatever the current directory happens to be.
    if nargin < 3 || isempty(outDir)
        outDir = fileparts(mfilename('fullpath'));
    end
    masterSeed = 42;
    rng(masterSeed, 'twister');
    
    bounds = readstruct(jsonFile);
    paramNames = fieldnames(bounds);
    dim = length(paramNames);
    
    unitSamples = lhsdesign(n, dim, 'criterion', 'maximin', 'iterations', 50);
    scenarios(1:n) = struct();
    
    for j = 1:dim
         name = paramNames{j};
         paramData = bounds.(name);
         
         if isfield(paramData, "type")
             distType = paramData.type;
         else
             distType = "uniform"; 
         end
         require_fields(name, paramData, distType);
         
        switch char(distType)
            case 'uniform'
                lower = paramData.lower;
                upper = paramData.upper;
                scaled = lower + (upper - lower).*unitSamples(:,j);
                
            case 'normal'
                mu = paramData.mean;
                sigma = paramData.std;
                scaled = norminv(unitSamples(:,j), mu, sigma);
                
            case 'lognormal'
                mu = paramData.mean;
                sigma = paramData.std;
                pd = makedist('Lognormal', 'mu', mu, 'sigma', sigma);
                scaled = icdf(pd, unitSamples(:,j));
                
            case 'weibull'
                a = paramData.scale;
                b = paramData.shape;
                pd = makedist('Weibull', 'a', a, 'b', b);
                scaled = icdf(pd, unitSamples(:,j));
                
            case 'beta'
                alpha = paramData.alpha;
                betaParam = paramData.beta;
                pd = makedist('Beta', 'a', alpha, 'b', betaParam);
                rawSamples = icdf(pd, unitSamples(:,j));
                lower = paramData.lower;
                upper = paramData.upper;
                scaled = lower + (upper - lower) .* rawSamples;
                
            case 'discrete'
                values = paramData.values;
                k = length(values);
                idx = floor(unitSamples(:,j)*k) + 1;
                idx(idx > k) = k;
                scaled = values(idx);
                
            otherwise
                error("Distribution type for %s not recognized.", name);
        end
        
        % Optional hard bounds clipping if specified (skip for beta since it scales directly to bounds)
        if ~strcmpi(string(distType), 'beta') && isfield(paramData, "lower") && isfield(paramData, "upper")
            scaled = max(scaled, paramData.lower);
            scaled = min(scaled, paramData.upper);
        end
        
        for i = 1:n
            scenarios(i).(name) = scaled(i);
        end
    end
    
    for i = 1:n
        scenarios(i).RunID = i;
        scenarios(i).Seed = masterSeed + i;
    end
    
    mcTable = struct2table(scenarios);
    writetable(mcTable, fullfile(outDir, 'mc_scenarios.csv'));
end

function require_fields(name, paramData, distType)
%REQUIRE_FIELDS Check one parameter's entry before sampling from it.
%   Catches a malformed bounds file here, naming the parameter, instead of
%   letting it surface later as a bare "reference to non-existent field".
    switch char(distType)
        case 'uniform',   needed = {'lower', 'upper'};
        case 'normal',    needed = {'mean', 'std'};
        case 'lognormal', needed = {'mean', 'std'};
        case 'weibull',   needed = {'scale', 'shape'};
        case 'beta',      needed = {'alpha', 'beta', 'lower', 'upper'};
        case 'discrete',  needed = {'values'};
        otherwise
            error('generateScenarios_ctv:badType', ...
                ['Parameter "%s" has type "%s". Use uniform, normal, ', ...
                 'lognormal, weibull, beta or discrete.'], name, char(distType));
    end

    missing = needed(~isfield(paramData, needed));
    if ~isempty(missing)
        error('generateScenarios_ctv:missingField', ...
            'Parameter "%s" (type %s) is missing: %s. It needs: %s.', ...
            name, char(distType), strjoin(missing, ', '), strjoin(needed, ', '));
    end
    if ~isfield(paramData, 'nominal')
        error('generateScenarios_ctv:missingNominal', ...
            'Parameter "%s" has no "nominal" value, which the nominal run uses.', name);
    end
    if all(isfield(paramData, {'lower', 'upper'})) && paramData.lower > paramData.upper
        error('generateScenarios_ctv:badBounds', ...
            'Parameter "%s" has lower (%g) above upper (%g).', ...
            name, paramData.lower, paramData.upper);
    end
end
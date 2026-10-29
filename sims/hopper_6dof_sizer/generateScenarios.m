function [scenarios, mcTable] = generateScenarios(jsonFile, n)
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
    writetable(mcTable, 'mc_scenarios.csv');
end
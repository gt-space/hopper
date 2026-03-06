function scenarios = generateScenarios(jsonFile, n)

    % Read JSON
    bounds = readstruct(jsonFile);
    paramNames = fieldnames(bounds);
    dim = length(paramNames);

    % Generate Latin Hypercube samples in [0,1] for Continous Parameters
    unitSamples = lhsdesign(n, dim, 'criterion', 'maximin', 'iterations', 50);

    % Preallocate struct array
    scenarios(1:n) = struct();

    for j = 1:dim
         name = paramNames{j};
         paramData = bounds.(name);

        % --- Continuous parameter ---
        if isfield(paramData, "lower")

            lower = paramData.lower;
            upper = paramData.upper;

            scaled = lower + (upper - lower).*unitSamples(:,j);

        % --- Discrete parameter ---
        elseif isfield(paramData, "values")

            values = paramData.values;
            k = length(values);

            % Map [0,1] to integer bins
            idx = floor(unitSamples(:,j)*k) + 1;

            % Fix edge case if sample == 1
            idx(idx > k) = k;

            scaled = values(idx);

        else
            error("Parameter %s not properly defined.", name);
        end

        for i = 1:n
            scenarios(i).(name) = scaled(i);
        end
    end
end
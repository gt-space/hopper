classdef VariantConfigTest < matlab.unittest.TestCase
    % variantConfig: defaults, choices by name, and typos rejected.

    methods (Test)
        function defaultsAreFirstChoices(tc)
            V = variantConfig();
            for f = fieldnames(variantAxes())'
                tc.verifyEqual(V.(['VARIANT_' f{1}]), 1, f{1});
            end
        end

        function choicesByName(tc)
            V = variantConfig('SLOSH', 'off', 'propulsion', "fluid_network");
            tc.verifyEqual(V.VARIANT_SLOSH, 2);
            tc.verifyEqual(V.VARIANT_PROPULSION, 2);
            tc.verifyEqual(V.VARIANT_GROUND, 1);
            tc.verifyEqual(variantConfig(struct('GROUND', 'strut')).VARIANT_GROUND, 2);
        end

        function typosAreErrors(tc)
            tc.verifyError(@() variantConfig('SLOSHH', 'off'), 'variantConfig:axis');
            tc.verifyError(@() variantConfig('SLOSH', 'none'), 'variantConfig:choice');
            tc.verifyError(@() variantConfig('SLOSH'), 'variantConfig:args');
        end
    end
end

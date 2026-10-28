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

        function ctvDefaults(tc)
            % the CTV burns no propellant: slosh off, constant mass properties
            V = variantConfig('VEHICLE', 'ctv');
            tc.verifyEqual(V.VARIANT_VEHICLE, 2);
            tc.verifyEqual(V.VARIANT_SLOSH, 2);
            tc.verifyEqual(V.VARIANT_MASS_PROPS, 2);
            V = variantConfig('VEHICLE', 'ctv', 'SLOSH', 'full');   % named axes win
            tc.verifyEqual(V.VARIANT_SLOSH, 1);
        end

        function ctvInputs(tc)
            % CAD mass and thrust; CAD Y (long axis) becomes body x
            old = hopperVehicle('ctv');
            tc.addTeardown(@() hopperVehicle(old));
            IN = mission_inputs();
            tc.verifyEqual(IN.mass.dry_override + IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass, 2.61488, 'AbsTol', 1e-9);
            tc.verifyEqual(IN.propulsion.nominal_thrust, 7.875 * 4.44822, 'AbsTol', 1e-9);
            tc.verifyEqual(IN.mass.dry_cg, 0.34921, 'AbsTol', 1e-9);
            tc.verifyEqual(IN.mass.dry_I(1, 1), 12321923.43e-9, 'RelTol', 1e-12);
            tc.verifyEqual(hopperVehicle(), 'ctv');
            hopperVehicle(old);
            tc.verifyEqual(mission_inputs().mass.dry_override, 116.86);
        end

        function typosAreErrors(tc)
            tc.verifyError(@() variantConfig('SLOSHH', 'off'), 'variantConfig:axis');
            tc.verifyError(@() variantConfig('SLOSH', 'none'), 'variantConfig:choice');
            tc.verifyError(@() variantConfig('SLOSH'), 'variantConfig:args');
        end
    end
end

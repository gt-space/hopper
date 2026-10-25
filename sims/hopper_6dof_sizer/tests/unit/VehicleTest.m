classdef VehicleTest < matlab.unittest.TestCase
    % Sizing, mass properties, slosh damping and attitude helpers.

    properties
        IN
        TANKS
    end

    methods (TestClassSetup)
        function size(tc)
            tc.IN = mission_inputs();
            tc.TANKS = size_tanks(tc.IN);
        end
    end

    methods (Test)
        function massPropertiesArePhysical(tc)
            MP = mass_properties_setup(tc.IN, tc.TANKS, 116.86);
            tc.verifyTrue(all(isfinite(MP.cg_init(:))));
            I = MP.MoI_init;
            tc.verifyEqual(I, I', 'AbsTol', 1e-9, 'inertia must be symmetric');
            tc.verifyGreaterThan(eig((I + I') / 2), 0, 'inertia must be positive definite');
            % CG sits between the engine mount and the vehicle top
            tc.verifyGreaterThan(MP.cg_init(1), tc.IN.mount_cg);
            tc.verifyLessThan(MP.cg_init(1), tc.IN.mission.max_vehicle_height);
        end

        function dispersionsMoveMassProperties(tc)
            nom = mass_properties_setup(tc.IN, tc.TANKS, 116.86);
            IN2 = tc.IN; IN2.inertia_factor = 1.05;
            hi = mass_properties_setup(IN2, tc.TANKS, 116.86);
            tc.verifyGreaterThan(hi.MoI_init(2, 2), nom.MoI_init(2, 2));
        end

        function dampingIsPositiveAndFinite(tc)
            B = load('baffle_lut.mat', 'ox_baffle', 'baffle_opts', 'tank_r', 'ox_density', 'ox_nu');
            m = linspace(0.1, tc.IN.propulsion.oxidizer_mass, 50);
            ob = B.ox_baffle;
            zl = lateralDampingCalc(m, B.ox_density, B.ox_nu, B.tank_r, ob.w, ob.Nb, ob.span, ob.t, B.baffle_opts.slosh_amp);
            za = axialDampingCalc(m, B.ox_density, B.ox_nu, B.tank_r, ob.w, ob.Nb, ob.span, ob.t, B.baffle_opts.slosh_amp);
            tc.verifyTrue(all(isfinite(zl)) && all(zl > 0));
            tc.verifyTrue(all(isfinite(za)) && all(za > 0));
        end

        function quaternionToEulerIdentity(tc)
            [phi, theta, psi] = Quaternion_to_Euler(1, 0, 0, 0);
            tc.verifyEqual([phi theta psi], [0 0 0], 'AbsTol', 1e-12);
        end

        function hopperRootFindsProject(tc)
            tc.verifyTrue(isfile(hopperRoot('HopperSim.prj')));
            tc.verifyTrue(isfolder(hopperRoot('data')));
        end
    end
end

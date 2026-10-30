classdef VehicleConfigTest < matlab.unittest.TestCase
    % Vehicle config files (vehicles/<vehicle>/<vehicle>_config.m): they
    % load and pass vehicleConfig's checks, the derived values are right,
    % the setup reads them, a geometry change reaches the mass properties,
    % and writeVehicleConfig round-trips a config exactly.

    properties (TestParameter)
        vehicle = cellstr(variantAxes().VEHICLE)   % every vehicle the VEHICLE axis lists
    end

    properties
        Dir   % temporary folder on the path for scratch config files
    end

    methods (TestMethodSetup)
        function scratchFolder(tc)
            f = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(f.Folder));
            tc.Dir = f.Folder;
        end
    end

    methods (Test)
        function configLoadsAndChecks(tc, vehicle)
            VC = vehicleConfig(vehicle);
            tc.verifyEqual(VC.name, vehicle);
            tc.verifyGreaterThan(VC.mass.wet, VC.mass.dry);
        end

        function derivedValues(tc)
            VC = vehicleConfig('hopper');
            tc.verifyEqual(VC.tanks.fu.floor_x, 1.0668 - 18.91*0.0254/2);
            tc.verifyEqual(VC.tanks.ox.floor_x, 1.854 - 29.18*0.0254/2);
            tc.verifyEqual(VC.mass.wet, 114.894 + 16 + 13);
            % CAD dry block: CAD Y from the footpad bottoms becomes body x; roll is CAD Lyy
            tc.verifyEqual(VC.mass.dry_cg, -0.169015 + 1.208957, 'AbsTol', 1e-12);
            tc.verifyEqual(diag(VC.mass.dry_I)', [12.6948 52.3340 51.5454]);
            tc.verifyEqual(VC.engine.tvc.natural_freq, 36, 'AbsTol', 1e-12);
            tc.verifyEqual(VC.structure.modal.k_engine, (0.3556 / 2.44)^2);
            d = VC.legs.foot_radius;   % legs 1-4 at +y -y +z -z, as the models expect
            tc.verifyEqual(VC.legs.r_foot, [0 d 0; 0 -d 0; 0 0 d; 0 0 -d]');
        end

        function setupReadsTheConfig(tc, vehicle)
            old = hopperVehicle(vehicle);
            restore = onCleanup(@() hopperVehicle(old));
            IN = mission_inputs();
            VC = vehicleConfig(vehicle);
            tc.verifyEqual(IN.vehicle.name, vehicle);
            tc.verifyEqual(IN.engine_cg, VC.engine.cg_x);
            tc.verifyEqual(IN.tanks.max_radius, VC.tanks.radius);
            tc.verifyEqual(IN.propulsion.nominal_thrust, VC.engine.nominal_thrust);
            V = vehicleSetup(IN);
            tc.verifyEqual(V.STRUCT.legs.k, VC.legs.stiffness);
            tc.verifyEqual(V.VEH.mass.dry, VC.mass.dry);
            tc.verifyEqual(V.VEH.k_engine, VC.structure.modal.k_engine);
        end

        function tankGeometryReachesMassProperties(tc)
            % a longer fuel tank (same centre) changes the propellant column,
            % so the full-tank CG and the CG tables move
            raw = hopper_config();
            raw.tanks.fu.height = raw.tanks.fu.height + 0.1;
            raw.name = 'tallfu';
            writeVehicleConfig(raw, fullfile(tc.Dir, 'tallfu_config.m'));
            IN0 = mission_inputs();
            IN1 = vehicleInputs(IN0, vehicleConfig('tallfu'));
            [mp0, V0] = VehicleConfigTest.massProps(IN0);
            [mp1, V1] = VehicleConfigTest.massProps(IN1);
            tc.verifyNotEqual(mp1.cg_init, mp0.cg_init, 'full-tank CG');
            tc.verifyNotEqual(mp1.dcg_dMfu_table, mp0.dcg_dMfu_table, 'fuel CG table');
            tc.verifyEqual(IN1.vehicle.tanks.fu.floor_x, IN0.vehicle.tanks.fu.floor_x - 0.05, 'AbsTol', 1e-12, ...
                'tank bottom the Slosh and Mass_Properties blocks read');
            tc.verifyEqual(V1.TANKS.singular.fuel.h, V0.TANKS.singular.fuel.h + 0.1, 'AbsTol', 1e-12);
        end

        function checksCatchBadGeometry(tc)
            raw = hopper_config();
            raw.tanks.ox.cg_x = raw.tanks.fu.cg_x;   % tanks on top of each other
            writeVehicleConfig(raw, fullfile(tc.Dir, 'overlap_config.m'));
            tc.verifyError(@() vehicleConfig('overlap'), 'vehicleConfig:geometry');
            raw = hopper_config();
            raw.height = 1;                          % the ox tank sticks out the top
            writeVehicleConfig(raw, fullfile(tc.Dir, 'short_config.m'));
            tc.verifyError(@() vehicleConfig('short'), 'vehicleConfig:geometry');
            raw = rmfield(hopper_config(), 'aero');
            writeVehicleConfig(raw, fullfile(tc.Dir, 'noaero_config.m'));
            tc.verifyError(@() vehicleConfig('noaero'), 'vehicleConfig:field');
            tc.verifyError(@() vehicleConfig('nosuchvehicle'), 'vehicleConfig:missing');
        end

        function writeRoundTrips(tc, vehicle)
            raw = feval([vehicle '_config']);
            writeVehicleConfig(raw, fullfile(tc.Dir, 'roundtrip_config.m'));
            tc.verifyEqual(roundtrip_config(), raw);
        end

        function sizerDraftsAConfig(tc)
            file = fullfile(tc.Dir, 'draft_config.m');
            [VC, S] = sizeVehicle('hopper', 'Write', file);
            tc.verifyTrue(isfile(file));
            tc.verifyEqual(VC.legs.stiffness, S.STRUCT.legs.k);
            tc.verifyEqual(VC.tanks.ox.mass, S.TANKS.singular.oxidizer.mass);
            tc.verifyEqual(draft_config(), VC, 'draft reads back as written');
            vehicleConfig('draft');   % and passes the checks
        end
    end

    methods (Static)
        function [MP, V] = massProps(IN)
            V = vehicleSetup(IN);
            MP = mass_properties_setup(IN, V.TANKS, V.OUT.Vehicle.DryMass);
        end
    end
end

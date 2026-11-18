classdef LandingModelTest < matlab.unittest.TestCase
    %LANDINGMODELTEST  Landing-loads model (src/analysis/landing).
    %   Geometry against the CAD, the nominal drop against structures'
    %   Ansys Flat11 case, crush-core and pad force caps, energy and step
    %   size. Runs cover only the first impact (tEnd 0.04 s).

    properties
        VC
    end

    methods (TestClassSetup)
        function load(tc)
            tc.VC = vehicleConfig('hopper');
        end
    end

    methods (Test)
        function geometryMatchesCad(tc)
            G = landingGear(Vehicle=tc.VC);
            tc.verifyEqual(G.L0(:, 1).', [0.906 0.909 0.909], 'AbsTol', 0.005);   % ~35.7 in pin to pin
            tc.verifyEqual(G.L0, repmat(G.L0(:, 1), 1, 4), 'AbsTol', 1e-9);       % four identical legs
            tc.verifyEqual(G.tip_angle_deg, 30.3, 'AbsTol', 0.5);
            % static: main ~5.1x, aux ~-2.2x the vertical foot load (pinned)
            N = G.u0(:, :, 1) \ [-1; 0; 0];
            tc.verifyEqual(N.', [5.11 -2.22 -2.22], 'AbsTol', 0.05);
        end

        function attitudeSetsWhichFootLandsFirst(tc)
            % drift toward leg 1, 2 deg tilt: leading -> foot 1 first,
            % braking -> foot 3 (trailing) first. Only that foot touches in 5 ms.
            G = landingGear(Vehicle=tc.VC);
            S = landingSurface('rigid');
            for a = {'leading', 1; 'braking', 3}.'
                C = landingCase(SinkRate=1.5, Lateral=0.5, Heading='leg', Tilt=2, ...
                    Attitude=a{1}, Vehicle=tc.VC);
                R = landingSim(G, S, C, tEnd=0.005, Log=false);
                [~, first] = max(R.peak.Fn);
                tc.verifyEqual(first, a{2}, a{1});
                tc.verifyEqual(nnz(R.peak.Fn > 0), 1, a{1});
            end
        end

        function nominalDropMatchesAnsys(tc)
            % Ansys Flat11: 145 kg, all four feet, main 45.2 kN, aux -18.1 kN
            R = landingSim(landingGear(Vehicle=tc.VC), landingSurface('rigid', mu=0), ...
                landingCase(SinkRate=1.5, Vehicle=tc.VC), tEnd=0.04, Log=false);
            tc.verifyEqual(max(R.peak.Nc(1, :)), 45.2e3, 'RelTol', 0.10);
            tc.verifyEqual(-min(R.peak.Nt(2:3, :), [], 'all'), 18.1e3, 'RelTol', 0.10);
        end

        function mainCoreCapsMainStrut(tc)
            G = landingGear(Vehicle=tc.VC, Core=crushCore('main', Plateau=12e3));
            R = landingSim(G, landingSurface('concrete'), landingCase(SinkRate=3.048, Vehicle=tc.VC), ...
                tEnd=0.04, Log=false);
            tc.verifyLessThanOrEqual(max(R.peak.Nc(1, :)), 12e3 * 1.001);
            tc.verifyGreaterThan(R.peak.stroke_core, 1e-3);
            tc.verifyGreaterThan(R.peak.E_core, 0);
        end

        function padCoreCapsGroundForce(tc)
            G = landingGear(Vehicle=tc.VC, Core=crushCore('pad', Plateau=4e3));
            R = landingSim(G, landingSurface('concrete'), landingCase(SinkRate=3.048, Vehicle=tc.VC), ...
                tEnd=0.04, Log=false);
            tc.verifyLessThanOrEqual(max(R.peak.Fn), 4e3 * 1.001);
            tc.verifyGreaterThan(R.peak.stroke_pad, 0.01);
        end

        function noEnergyCreated(tc)
            R = landingSim(landingGear(Vehicle=tc.VC), landingSurface('rigid', mu=0), ...
                landingCase(SinkRate=1.5, Vehicle=tc.VC), tEnd=0.04, Log=false);
            % elastic legs store energy at the end; it can't exceed the input
            tc.verifyLessThanOrEqual(R.final.E, R.peak.E_in * 1.001);
        end

        function softSoilDentsAndCushions(tc)
            G = landingGear(Vehicle=tc.VC);
            C = landingCase(SinkRate=1.5, Vehicle=tc.VC);
            Rs = landingSim(G, landingSurface('soil_soft'), C, tEnd=0.04, Log=false);
            Rc = landingSim(G, landingSurface('concrete'), C, tEnd=0.04, Log=false);
            tc.verifyGreaterThan(Rs.peak.dent, 0);
            tc.verifyLessThan(max(Rs.peak.Nc(1, :)), max(Rc.peak.Nc(1, :)));
        end

        function stepSizeConverged(tc)
            G = landingGear(Vehicle=tc.VC, Core=crushCore('main', Plateau=15e3));
            C = landingCase(SinkRate=3.048, Lateral=0.6, Tilt=2, Vehicle=tc.VC);
            S = landingSurface('concrete');
            a = landingSim(G, S, C, dt=2e-5, tEnd=0.04, Log=false);
            b = landingSim(G, S, C, dt=1e-5, tEnd=0.04, Log=false);
            tc.verifyEqual(a.peak.stroke_core, b.peak.stroke_core, 'RelTol', 0.03);
            tc.verifyEqual(min(a.peak.Nt, [], 'all'), min(b.peak.Nt, [], 'all'), 'RelTol', 0.03);
        end
    end
end

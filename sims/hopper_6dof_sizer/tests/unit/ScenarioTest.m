classdef ScenarioTest < matlab.unittest.TestCase
    % Monte Carlo scenario generation and the mission pass/fail check.

    methods (Test)
        function scenariosAreSeeded(tc)
            a = generateScenarios('mc_params.json', 10);
            b = generateScenarios('mc_params.json', 10);
            tc.verifyEqual(a, b);
            tc.verifyNumElements(a, 10);
        end

        function scenariosStayInBounds(tc)
            P = readstruct('mc_params.json');
            s = generateScenarios('mc_params.json', 25);
            for f = fieldnames(P)'
                v = [s.(f{1})];
                tc.verifyGreaterThanOrEqual(v, P.(f{1}).lower, f{1});
                tc.verifyLessThanOrEqual(v, P.(f{1}).upper, f{1});
            end
        end

        function goodFlightPasses(tc)
            r = checkConstraints(fakeResult());
            tc.verifyTrue(r.status.pass);
            tc.verifyEmpty(r.status.failed);
        end

        function eachLimitIsChecked(tc)
            r = fakeResult(); r.flight.landing_vel = 2;
            tc.verifyEqual(checkConstraints(r).status.failed, {'landing_vel'});
            r = fakeResult(); r.timeseries.x(end) = 3;
            tc.verifyEqual(checkConstraints(r).status.failed, {'landing_offset'});
            r = fakeResult(); r.timeseries.pitch(5) = 70;
            tc.verifyEqual(checkConstraints(r).status.failed, {'max_tilt_deg'});
            r = fakeResult(); r.flight.max_altitude = 40;
            tc.verifyEqual(checkConstraints(r).status.failed, {'max_altitude'});
        end

        function limitsCanBeOverridden(tc)
            r = fakeResult(); r.flight.landing_vel = 2;
            r = checkConstraints(r, struct('max_landing_vel', 3));
            tc.verifyTrue(r.status.pass);
        end
    end
end

function r = fakeResult()
r.vehicle = struct('twr_initial', 1.5, 'twr_final', 2, 'dry_mass', 117, 'wet_mass', 146);
r.flight = struct('max_altitude', 52, 'max_ascent_vel', 5, 'max_descent_vel', 5, 'landing_vel', 0.3);
r.timeseries = struct('x', zeros(10, 1), 'y', zeros(10, 1), 'pitch', 90 * ones(10, 1));
r.status = struct('success', true, 'error', '');
end

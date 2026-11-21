classdef ControllerWarmupTest < matlab.unittest.TestCase
    % The GMEKF warm-up delays references without modifying controller logic.

    methods (Test)
        function delaysScheduleAndTrajectorySignals(tc)
            W.tgrid = [0; 0.5; 1];
            W.unom = [9.8 0 0 0; 10 0 0 0; 9.8 0 0 0];
            W.T_ts = timeseries([9.8; 10; 9.8], W.tgrid);
            W.z_ref_ts = timeseries([0; 1; 0], W.tgrid);
            W.scenario = Simulink.SimulationData.Dataset;
            W.scenario = W.scenario.addElement(W.T_ts, 'thrust_ref');

            W = applyControllerWarmup(W, 10);

            tc.verifyEqual(W.controller_start_time, 10);
            tc.verifyEqual(W.tgrid, [10; 10.5; 11]);
            tc.verifyEqual(W.T_ts.Time, [10; 10.5; 11]);
            tc.verifyEqual(W.z_ref_ts.Time, [10; 10.5; 11]);
            tc.verifyEqual(W.scenario.getElement(1).Time, [10; 10.5; 11]);
            tc.verifyEqual(W.u_hold_command, zeros(4, 1));
        end

        function changingWarmupRestoresOriginalTimes(tc)
            W.tgrid = [0; 0.5; 1];
            W.unom = ones(3, 4);
            W.T_ts = timeseries([1; 1; 1], W.tgrid);

            W = applyControllerWarmup(W, 10);
            W = applyControllerWarmup(W, 0);

            tc.verifyEqual(W.tgrid, [0; 0.5; 1]);
            tc.verifyEqual(W.T_ts.Time, [0; 0.5; 1]);
            tc.verifyEqual(W.controller_start_time, 0);
        end

        function rejectsInvalidWarmup(tc)
            W.tgrid = [0; 1];
            W.unom = ones(2, 4);
            tc.verifyError(@() applyControllerWarmup(W, -1), 'applyControllerWarmup:invalidTime');
        end
    end
end

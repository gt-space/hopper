classdef FdirFlightTest < matlab.unittest.TestCase
    % Full flights of models/hopper_6dof.slx with FDIR and fault injection
    % on: an injected fault is detected and gets the GNC PDR response.
    %   lidar beam loss  -> R0, the beam is masked and the flight is nominal
    %   gyro stuck       -> R4 is not allowed above h_cut, so the vehicle
    %                       flies an abort to land, then cuts off at h_cut

    properties (Constant)
        Model = 'hopper_6dof'
        V = {'FDIR', 'on', 'FAULT_INJECTION', 'on'}
    end

    properties
        W
    end

    methods (TestClassSetup)
        function setup(tc)
            evalin('base', 'clearvars');
            tc.W = hopperWorkspace();
        end
    end

    methods (Test)
        function lidarBeamLossIsMasked(tc)
            out = tc.fly({'lidar', 'loss', 20, inf, 0, 2});
            [~, info] = fdir_load_catalog();
            st = FdirFlightTest.logged(out, 'fdir_status');
            beam2 = find(info.name == "Lidar beam 2");
            tc.verifyTrue(any(st(:, beam2) == 2), 'beam 2 loss never latched');
            tc.verifyEqual(max(FdirFlightTest.logged(out, 'fdir_level')), 0, 'response beyond R0');   % -1 none, 0 mask
            M = fdir_modes();
            mode = FdirFlightTest.logged(out, 'fdir_mode');
            % on the ground: the sim stops at contact, which can come a frame
            % or two before touchdown confirmation reaches LANDED
            tc.verifyTrue(ismember(double(mode(end)), double([M.TOUCHDOWN_CONFIRM M.LANDED])), 'final mode');
        end

        function gyroStuckAbortsToLand(tc)
            out = tc.fly({'gyro', 'stuck', 8});
            M = fdir_modes();
            mode = FdirFlightTest.logged(out, 'fdir_mode');
            tc.verifyTrue(any(mode == double(M.ABORT_LAND)), 'no abort to land');
            z = squeeze(out.get('z').Data);
            % the latched R4 may only cut off below h_cut
            tm = out.get('fdir_mode').Time;
            i_safe = find(mode == double(M.SAFED), 1);
            if ~isempty(i_safe)
                h = -interp1(out.get('z').Time, z, tm(i_safe));
                tc.verifyLessThan(h, fdir_mode_params(tc.W.IN).h_cut + 0.1, 'cut off above h_cut');
            end
            v = squeeze(out.get('velocity').Data);
            tc.verifyLessThan(-z(end), 0.1, 'did not come down');
            % R4 is latched, so the engine cuts at h_cut (PDR) and the
            % last 1.5 m are unpowered: ~3 m/s at touchdown
            tc.verifyLessThan(abs(v(end, 3)), 3.5, 'touchdown sink rate (m/s)');
        end
    end

    methods
        function out = fly(tc, varargin)
            W = tc.W;
            W.faults = faultScenario(varargin{:});
            out = sim(buildSimInput(tc.Model, [], W, variantConfig(tc.V{:})));
        end
    end

    methods (Static)
        function x = logged(out, name)
            % logged signal as time x components, double
            s = out.get(name);
            x = double(squeeze(s.Data));
            if size(x, 1) ~= numel(s.Time), x = x'; end
        end
    end
end

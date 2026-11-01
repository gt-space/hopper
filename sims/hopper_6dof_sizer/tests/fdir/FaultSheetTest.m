classdef FaultSheetTest < matlab.unittest.TestCase
    % inputs/fdir_fault_conditions.csv (FaultConditions.xlsx, from #31)
    % and the catalog's sheet_id links to it. Where the catalog's response
    % and the sheet's tier disagree is reported by fdir_sheet_check, not
    % failed here: that is the FDIR team's call (docs/FDIR_ARCHITECTURE.md).

    methods (Test)
        function sheetIsTheFullFaultList(tc)
            S = readtable(hopperRoot('inputs', 'fdir_fault_conditions.csv'), 'TextType', 'string', ...
                'VariableNamingRule', 'preserve');
            tc.verifyEqual(height(S), 102);
            tc.verifyEqual(numel(unique(S.ID)), 102, 'sheet IDs are unique');
            tc.verifyEqual(S.ID(1), "F-01");
        end

        function everyLinkResolves(tc)
            % fdir_sheet_check errors on a sheet_id the sheet doesn't have
            T = fdir_sheet_check();
            [~, C] = fdir_load_catalog();
            tc.verifyEqual(height(T), height(C));
            tc.verifyTrue(all(ismember(T.verdict, ["agrees", "differs", "conditional", "no sheet row"])));
        end

        function verdicts(tc)
            % R<k> against the sheet's tier; conditional responses have none
            T = fdir_sheet_check();
            v = @(id) T.verdict(T.id == id);
            tc.verifyEqual(v(12), "agrees", 'lidar all beams: R2, sheet tier 2');
            tc.verifyEqual(v(1), "differs", 'IMU loss: R4, sheet tier 2');
            tc.verifyEqual(v(19), "conditional", 'attitude envelope: R4 above h_cut');
            tc.verifyEqual(v(33), "no sheet row");
        end
    end
end

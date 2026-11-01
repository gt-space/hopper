function T = fdir_sheet_check(catalogFile, sheetFile)
% FDIR_SHEET_CHECK  Catalog responses against the team's fault sheet.
%
%   T = fdir_sheet_check()   one row per catalog fault
%
%   inputs/fdir_fault_conditions.csv is the "Fault List" sheet of
%   FaultConditions.xlsx (102 failure conditions over every subsystem,
%   with severity and ResponseTier = severity - 1), as #31 exported it.
%   inputs/fdir_fault_catalog.csv says how each monitor detects and
%   responds; its sheet_id column names the sheet row it answers to.
%
%   T columns: id, name, response (catalog), sheet_id, sheet_tier,
%   condition (sheet text), verdict:
%     "agrees"       catalog R<k> equals the sheet's tier
%     "differs"      it doesn't: the FDIR team decides which is right
%     "conditional"  NO_ARM, R4_ABOVE_HCUT, STUCK_THRUST: no single tier
%     "no sheet row" the catalog row has no sheet_id
%
%   Nothing here changes a response; it lists where the two disagree.

if nargin < 1, catalogFile = hopperRoot('inputs', 'fdir_fault_catalog.csv'); end
if nargin < 2, sheetFile = hopperRoot('inputs', 'fdir_fault_conditions.csv'); end
[~, C] = fdir_load_catalog(catalogFile);
S = readtable(sheetFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');

n = height(C);
tier = nan(n, 1); cond = strings(n, 1); verdict = strings(n, 1);
for k = 1:n
    sid = C.sheet_id(k);
    if sid == ""
        verdict(k) = "no sheet row";
        continue
    end
    j = find(S.ID == sid, 1);
    if isempty(j)
        error('fdir_sheet_check:id', 'Catalog fault %d names %s, which the sheet has no row for.', C.id(k), sid);
    end
    tier(k) = str2double(string(S.ResponseTier(j)));
    cond(k) = S.FailureCondition(j);
    r = regexp(C.response(k), '^R(\d)$', 'tokens', 'once');
    if isempty(r)
        verdict(k) = "conditional";
    elseif str2double(r{1}) == tier(k)
        verdict(k) = "agrees";
    else
        verdict(k) = "differs";
    end
end
T = table(C.id, C.name, C.response, C.sheet_id, tier, cond, verdict, ...
    'VariableNames', {'id', 'name', 'response', 'sheet_id', 'sheet_tier', 'condition', 'verdict'});
end

function make_fdir_catalog(csvPath, outPath)
% MAKE_FDIR_CATALOG  Regenerate fdir_fault_catalog_data.m from the CSV.
%
% fdir_fault_catalog.csv is the "Fault List" sheet of FaultConditions.xlsx
% in numeric F-ID order, plus ResponseTier (= severity - 1) and
% CourseOfAction (what that tier commands in fdir_response_mapping.m).
% The CSV is the source of truth for the catalog.
%
% A MATLAB Function block cannot read a CSV at runtime, so the numbers have
% to be baked into a .m file that code generation can fold to constants.
% That file is generated, never hand-edited:
%
%   make_fdir_catalog()                % rewrites fdir_fault_catalog_data.m
%   make_fdir_catalog(csv, out)        % explicit paths, used by the tests
%
% FDIR_Tests regenerates into a temporary file and compares, so a CSV edit
% that was never regenerated fails the suite.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(csvPath), csvPath = fullfile(here, 'fdir_fault_catalog.csv'); end
if nargin < 2 || isempty(outPath), outPath = fullfile(here, 'fdir_fault_catalog_data.m'); end

T = readtable(csvPath, 'TextType', 'string', 'VariableNamingRule', 'preserve');
n = height(T);

severity = zeros(1, n);
tier     = zeros(1, n);
for i = 1:n
    severity(i) = localNum(T.Severity(i));
    tier(i)     = localNum(T.ResponseTier(i));
end

% The sheet carries no course-of-action column; tier is always severity - 1,
% with an N/A severity (a fault that is not on the vehicle) pinned to R0.
expected = max(0, severity - 1);
if ~isequal(tier, expected)
    bad = find(tier ~= expected);
    error('make_fdir_catalog:TierMismatch', ...
        'ResponseTier is not severity - 1 for %s', strjoin("F-" + string(bad), ', '));
end

lines = strings(0, 1);
lines(end+1) = "function [severity, tier] = fdir_fault_catalog_data()";
lines(end+1) = "% FDIR_FAULT_CATALOG_DATA  Severity and response tier for all " + n + " faults.";
lines(end+1) = "%";
lines(end+1) = "% GENERATED FILE - do not edit by hand.";
lines(end+1) = "% Run make_fdir_catalog to rebuild it from fdir_fault_catalog.csv,";
lines(end+1) = "% which is the ""Fault List"" sheet of FaultConditions.xlsx.";
lines(end+1) = "%";
lines(end+1) = "% Severity 0 means the sheet lists the fault as N/A (not on the vehicle).";
lines(end+1) = "% Tier is severity - 1: 0 -> R0, 1 -> R1, 2 -> R2, 3 -> R3, 4 -> R4.";
lines(end+1) = "%";
lines(end+1) = "%#codegen";
lines(end+1) = "";
lines(end+1) = "severity = [ ...";
for i = 1:n
    lines(end+1) = sprintf('    %d, ...   %s', severity(i), localComment(T, i, tier(i))); %#ok<AGROW>
end
lines(end+1) = "    ];";
lines(end+1) = "";
lines(end+1) = "tier = [ ...";
for r = 1:10:n
    last = min(r + 9, n);
    nums = strjoin(compose("%d", tier(r:last)), ", ");
    if last < n, nums = nums + ","; end
    lines(end+1) = sprintf('    %s ...   %% %s..%s', nums, localId(T, r), localId(T, last)); %#ok<AGROW>
end
lines(end+1) = "    ];";
lines(end+1) = "end";
lines(end+1) = "";

fid = fopen(outPath, 'w');
if fid < 0
    error('make_fdir_catalog:Open', 'Could not open %s for writing.', outPath);
end
closer = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines(i));
end

fprintf('Wrote %s from %s (%d faults).\n', outPath, csvPath, n);
end

function v = localNum(x)
% The sheet writes N/A for a fault that is not on the vehicle.
s = strtrim(string(x));
if s == "" || ismissing(s) || upper(s) == "N/A" || upper(s) == "NA"
    v = 0;
else
    v = double(str2double(s));
end
end

function id = localId(T, i)
id = strtrim(string(T.ID(i)));
end

function c = localComment(T, i, tierValue)
% e.g. % F-15 GNC / Control: Loss of control authority (TWR < 1) -> R4
cond = localAscii(strtrim(string(T.FailureCondition(i))));
sub  = localAscii(strtrim(string(T.Subsystem(i))));
c = char(sprintf('%% %s %s: %s -> R%d', localId(T, i), sub, cond, tierValue));
end

function s = localAscii(s)
% Keep the generated file plain ASCII so it regenerates byte-identically
% whatever the locale. The CSV keeps the exact sheet text.
s = replace(s, char(177), "+/-");    % plus-minus
s = replace(s, char(8211), "-");     % en dash
s = replace(s, char(8212), "-");     % em dash
s = replace(s, char(8217), "'");     % right single quote
s = replace(s, char(8220), """");    % left double quote
s = replace(s, char(8221), """");    % right double quote
s = regexprep(s, '[^\x00-\x7F]', '');
end

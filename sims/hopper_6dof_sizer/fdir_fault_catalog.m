function T = fdir_fault_catalog(csvPath)
% FDIR_FAULT_CATALOG  The fault catalog as a readable table.
%
% For people, not for the model. The MATLAB Function blocks use
% fdir_fault_catalog_data (numeric, code-generatable); this reads the CSV
% so you can see what a fault ID actually means:
%
%   T = fdir_fault_catalog();
%   T(T.Severity == 5, :)                  % everything that commands R4
%   T(T.ID == "F-15", :).CourseOfAction
%
% Init_FDIR publishes it to the base workspace as fault_catalog, so the
% full table is in the workspace while the sim runs.
%
% Not codegen: it returns text.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(csvPath)
    csvPath = fullfile(here, 'fdir_fault_catalog.csv');
end

raw = readtable(csvPath, 'TextType', 'string', 'VariableNamingRule', 'preserve');

n = height(raw);
severity = zeros(n, 1);
tier     = zeros(n, 1);
for i = 1:n
    severity(i) = localNum(raw.Severity(i));
    tier(i)     = localNum(raw.ResponseTier(i));
end

T = table( ...
    strtrim(string(raw.ID)), ...
    strtrim(string(raw.Subsystem)), ...
    strtrim(string(raw.FailureCondition)), ...
    severity, ...
    tier, ...
    strtrim(string(raw.CourseOfAction)), ...
    'VariableNames', {'ID', 'Subsystem', 'FailureCondition', ...
                      'Severity', 'ResponseTier', 'CourseOfAction'});
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

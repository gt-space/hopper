function [cat, info] = fdir_load_catalog(file)
% FDIR_LOAD_CATALOG  Fault table (inputs/fdir_fault_catalog.csv) as numbers.
%
%   [cat, info] = fdir_load_catalog()
%
%   cat   N x 11 double, one row per fault, columns (fdir_catalog_cols):
%         n, m, cusum_limit, filter, clear_policy, phase_mask, response,
%         escalate_to, escalate_after_s, inhibited, target
%         so a MATLAB Function block can take it as a parameter.
%   info  the table as read (names, FTA rows, notes) for logs and reports.
%
% Format (FDIR CoDR slide 15, GNC PDR slide 78): id maps to the fault-tree
% row; latch when N of the last M evaluations fail (N_OF_M), on N in a row
% (CONSECUTIVE), or when the monitor's CUSUM passes cusum_limit (CUSUM);
% clear policy NEVER_IN_FLIGHT / AUTO_CLEAR / GROUND_ONLY; phase mask from
% PAD STARTUP ASCENT HOVER DESCENT LANDED (or ALL); response R0..R4, NO_ARM
% (IBIT), or the conditional R4_ABOVE_HCUT and STUCK_THRUST.

if nargin < 1, file = hopperRoot('inputs', 'fdir_fault_catalog.csv'); end
info = readtable(file, 'TextType', 'string', 'Delimiter', ',');
info = fillmissing(info, 'constant', "", 'DataVariables', @isstring);
K = fdir_catalog_cols();
M = fdir_modes();
n = height(info);
if ~isequal(info.id(:)', 1:n)
    error('fdir_load_catalog:ids', 'Fault ids must run 1..%d in order.', n);
end

cat = zeros(n, K.NCOL);
cat(:, K.N) = info.n;
cat(:, K.M) = info.m;
cat(:, K.CUSUM) = info.cusum_limit;
cat(:, K.FILTER) = code(info.filter, ["N_OF_M", "CONSECUTIVE", "CUSUM"], "filter");
cat(:, K.CLEAR) = code(info.clear_policy, ["NEVER_IN_FLIGHT", "AUTO_CLEAR", "GROUND_ONLY"], "clear_policy");
phases = ["PAD", "STARTUP", "ASCENT", "HOVER", "DESCENT", "LANDED"];
for k = 1:n
    p = strtrim(split(info.phase_mask(k), "|"));
    if any(p == "ALL")
        cat(k, K.PHASE_MASK) = double(M.ALL_PHASES);
    else
        bits = code(p, phases, "phase_mask");
        cat(k, K.PHASE_MASK) = sum(2 .^ (bits - 1));
    end
end
responses = ["R0", "R1", "R2", "R3", "R4", "NO_ARM", "R4_ABOVE_HCUT", "STUCK_THRUST"];
cat(:, K.RESPONSE) = code(info.response, responses, "response") - 1;   % R0 = 0 ... R4 = 4, NO_ARM 5, ...
esc = info.escalate_to;
has = strlength(esc) > 0;
cat(:, K.ESC_TO) = -1;
cat(has, K.ESC_TO) = code(esc(has), responses(1:5), "escalate_to") - 1;
cat(:, K.ESC_AFTER) = info.escalate_after_s;
cat(isnan(cat(:, K.ESC_AFTER)), K.ESC_AFTER) = inf;
cat(:, K.INHIBITED) = info.inhibited;
cat(:, K.TARGET) = code(info.target, fdir_targets(), "target");
end

function c = code(values, names, column)
values = string(values(:));
c = zeros(numel(values), 1);
for k = 1:numel(values)
    hit = find(names == values(k), 1);
    if isempty(hit)
        error('fdir_load_catalog:value', 'Unknown %s ''%s''. Allowed: %s.', column, values(k), strjoin(names, ', '));
    end
    c(k) = hit;
end
end

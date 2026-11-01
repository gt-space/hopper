function full_fails = fdir_expand_raw_fails(monitor_fails, fault_map, num_faults)
% FDIR_EXPAND_RAW_FAILS  Monitor checks -> catalog-indexed fault vector.
%
% The four sensor monitors produce 4 sensors x 4 checks = 16 flags, but
% persistence and fault_config are indexed by catalog fault ID over 102
% rows. fault_map (see fdir_config) gives the catalog ID of each check.
%
%#codegen

full_fails = false(num_faults, 1);
mf = logical(monitor_fails(:));

for k = 1:numel(mf)
    id = double(fault_map(k));
    full_fails(id) = full_fails(id) || mf(k);
end
end

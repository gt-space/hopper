function p_state = fdir_persistence_state(num_items, m_max)
% FDIR_PERSISTENCE_STATE  Fresh N-of-M filter state for num_items monitors.
%
% m_max is the widest M window any monitor uses; each one reads back only
% the last M_Window columns of its row.
%
%#codegen

p_state = struct( ...
    'history',          false(num_items, m_max), ...
    'states_vec',       zeros(num_items, 1, 'uint8'), ...
    'first_timestamps', zeros(num_items, 1), ...
    'sim_time',         0.0);
end

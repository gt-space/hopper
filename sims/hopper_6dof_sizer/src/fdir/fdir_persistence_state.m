function p_state = fdir_persistence_state(num_faults, m_max)
% FDIR_PERSISTENCE_STATE  Fresh N-of-M filter state for num_faults faults.
%
% m_max is the widest M window any fault uses; each fault reads back only
% the last M_Window columns of its row.
%
%#codegen

p_state = struct( ...
    'history',          false(num_faults, m_max), ...
    'states_vec',       zeros(num_faults, 1, 'uint8'), ...
    'first_timestamps', zeros(num_faults, 1), ...
    'sim_time',         0.0);
end

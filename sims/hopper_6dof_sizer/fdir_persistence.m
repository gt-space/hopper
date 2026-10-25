function [states, fail_counts, p_state] = fdir_persistence( ...
        raw_fails, ground_resets, flight_state, check_config, p_state)
% FDIR_PERSISTENCE  Per-check N-of-M persistence filter.
%
% N and M come from check_config(i).N_Threshold / .M_Window rather than a
% hard-coded 4-of-10, so loss and stuck checks can latch inside one frame.
%
% This runs over the 16 monitor checks, not the 102 catalog rows: checks
% that share a catalog row do not share a threshold (Barometer loss needs
% N = 1, Barometer drift needs N = 4, both report against F-04). Use
% fdir_fold_checks to turn the result into catalog-indexed states.
%
%   states       num_checks x 1 struct array matching Fault_State_Bus
%   fail_counts  num_checks x 1 uint8, failures inside each check's window
%   p_state      filter state, feed back in next frame
%                (see fdir_persistence_state)
%
%#codegen

ACTIVE_PHASES = [2, 3, 4];   % Ascent, Burn, Descent
INACTIVE = uint8(0);
PENDING  = uint8(1);
LATCHED  = uint8(2);

num_items = numel(raw_fails);
m_max = size(p_state.history, 2);

p_state.sim_time = p_state.sim_time + 0.01;   % 100 Hz frame rate increment

fail_counts = zeros(num_items, 1, 'uint8');

% Output bus struct template matching Fault_State_Bus
single_fault_struct = struct( ...
    'Fault_ID',              uint8(0), ...
    'Monitor_Valid',         true, ...
    'Raw_Fail',              false, ...
    'Persistence_State',     uint8(0), ...
    'Fail_Count',            uint8(0), ...
    'First_Fault_Timestamp', 0.0);
states = repmat(single_fault_struct, num_items, 1);

in_active_phase = any(flight_state == ACTIVE_PHASES);

for i = 1:num_items
    states(i).Fault_ID = uint8(i);
    states(i).Monitor_Valid = true;
    states(i).Raw_Fail = logical(raw_fails(i));

    n_thresh = double(check_config(i).N_Threshold);
    m_window = min(double(check_config(i).M_Window), m_max);

    % Handle individual ground reset
    if ground_resets(i)
        p_state.history(i, :) = false(1, m_max);
        p_state.states_vec(i) = INACTIVE;
        p_state.first_timestamps(i) = 0.0;
        fail_counts(i) = uint8(0);
        states(i).Persistence_State = p_state.states_vec(i);
        states(i).Fail_Count = fail_counts(i);
        states(i).First_Fault_Timestamp = p_state.first_timestamps(i);
        continue;
    end

    % Pause evaluation outside active phases
    if ~in_active_phase
        if p_state.states_vec(i) ~= LATCHED
            p_state.history(i, :) = false(1, m_max);
            p_state.states_vec(i) = INACTIVE;
        end
        fail_counts(i) = local_window_count(p_state.history, i, m_window);
        states(i).Persistence_State = p_state.states_vec(i);
        states(i).Fail_Count = fail_counts(i);
        states(i).First_Fault_Timestamp = p_state.first_timestamps(i);
        continue;
    end

    % Hold state if already latched
    if p_state.states_vec(i) == LATCHED
        fail_counts(i) = local_window_count(p_state.history, i, m_window);
        states(i).Persistence_State = p_state.states_vec(i);
        states(i).Fail_Count = fail_counts(i);
        states(i).First_Fault_Timestamp = p_state.first_timestamps(i);
        continue;
    end

    % Shift rolling window and insert sample
    p_state.history(i, :) = [p_state.history(i, 2:end), logical(raw_fails(i))];
    fail_counts(i) = local_window_count(p_state.history, i, m_window);

    % State transitions. A new failure enters PENDING and is then tested
    % against N in the same frame, so N = 1 latches on the first failing
    % frame and N = 4 still latches on the fourth.
    if p_state.states_vec(i) == INACTIVE && raw_fails(i)
        p_state.states_vec(i) = PENDING;
        p_state.first_timestamps(i) = p_state.sim_time;
    end

    if p_state.states_vec(i) == PENDING
        if double(fail_counts(i)) >= n_thresh
            p_state.states_vec(i) = LATCHED;
        elseif fail_counts(i) == 0
            p_state.states_vec(i) = INACTIVE;
            p_state.first_timestamps(i) = 0.0;
        end
    end

    states(i).Persistence_State = p_state.states_vec(i);
    states(i).Fail_Count = fail_counts(i);
    states(i).First_Fault_Timestamp = p_state.first_timestamps(i);
end
end

function c = local_window_count(history, i, m_window)
% Failures inside this fault's own M window, counting back from the newest
% sample at the right-hand end of the row. Written as a fixed-bound loop
% rather than a variable-length slice so it compiles in a MATLAB Function
% block without dynamic memory allocation.
%#codegen
m_max = size(history, 2);
first = m_max - m_window;
c = uint8(0);
for j = 1:m_max
    if j > first && history(i, j)
        c = c + uint8(1);
    end
end
end

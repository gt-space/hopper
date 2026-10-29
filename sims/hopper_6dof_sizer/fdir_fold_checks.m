function fault_states = fdir_fold_checks(check_states, fault_map, num_faults)
% FDIR_FOLD_CHECKS  Latched monitor checks -> catalog-indexed fault states.
%
% Persistence runs over the 16 monitor checks, because checks sharing a
% catalog row do not share a threshold. Arbitration runs over the 102
% catalog rows, because that is where Response_Tier lives. This folds one
% into the other: a catalog row carries the worst persistence state of any
% check mapped to it, so the row is latched if any of its checks latched.
%
% Rows with no monitor behind them stay INACTIVE.
%
%   check_states  num_checks x 1 struct array (Fault_State_Bus)
%   fault_states  num_faults x 1 struct array (Fault_State_Bus)
%
%#codegen

single_fault_struct = struct( ...
    'Fault_ID',              uint8(0), ...
    'Monitor_Valid',         true, ...
    'Raw_Fail',              false, ...
    'Persistence_State',     uint8(0), ...
    'Fail_Count',            uint8(0), ...
    'First_Fault_Timestamp', 0.0);
fault_states = repmat(single_fault_struct, num_faults, 1);

for i = 1:num_faults
    fault_states(i).Fault_ID = uint8(i);
end

for k = 1:numel(check_states)
    id = double(fault_map(k));

    % Worst state wins: INACTIVE(0) < PENDING(1) < LATCHED(2).
    if check_states(k).Persistence_State > fault_states(id).Persistence_State
        fault_states(id).Persistence_State = check_states(k).Persistence_State;
    end

    if check_states(k).Fail_Count > fault_states(id).Fail_Count
        fault_states(id).Fail_Count = check_states(k).Fail_Count;
    end

    fault_states(id).Raw_Fail = fault_states(id).Raw_Fail || check_states(k).Raw_Fail;
    fault_states(id).Monitor_Valid = fault_states(id).Monitor_Valid && check_states(k).Monitor_Valid;

    % Earliest non-zero timestamp across the row's checks.
    ts = check_states(k).First_Fault_Timestamp;
    if ts > 0 && (fault_states(id).First_Fault_Timestamp == 0 || ...
                  ts < fault_states(id).First_Fault_Timestamp)
        fault_states(id).First_Fault_Timestamp = ts;
    end
end
end

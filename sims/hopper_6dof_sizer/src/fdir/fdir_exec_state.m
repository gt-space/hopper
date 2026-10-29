function e = fdir_exec_state(n_faults)
% FDIR_EXEC_STATE  Initial state for fdir_executive.
%
%#codegen

H = 16;                                   % longest M window supported
e.history    = false(n_faults, H);        % last H pass/fail bits per fault
e.consec     = zeros(n_faults, 1);        % consecutive failing evaluations
e.cusum      = zeros(n_faults, 1);
e.status     = zeros(n_faults, 1, 'uint8');   % 0 INACTIVE, 1 PENDING, 2 LATCHED
e.fail_count = zeros(n_faults, 1);        % total failing evaluations this flight
e.first_t    = -ones(n_faults, 1);        % time first latched (s), -1 never
e.level      = int8(-1);                  % latched response level, -1 none (only escalates)
e.first_fault = 0;                        % id of the first fault to latch
e.first_fault_t = -1;
end

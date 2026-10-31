function st = fdir_frame_state(n_faults, n_fault_rows, dt)
% FDIR_FRAME_STATE  Initial state for fdir_frame.
%
%   n_faults      fault-table (catalog) rows
%   n_fault_rows  fault-injection table rows (faultScenario); unused since
%   dt            the sensors and their faults moved to sensor_frame, kept
%                 so callers don't change
%
%#codegen

st.mon   = gnc_monitor_state();
st.exec  = fdir_exec_state(n_faults);
st.modes = fdir_mode_state();
st.abort = fdir_abort_state();
st.prev_f = 9.80665;
end

function st = sensor_frame_state(n_fault_rows, dt)
%SENSOR_FRAME_STATE  Initial state for sensor_frame.
%   n_fault_rows  fault-injection table rows (faultScenario)
%   dt            the Sensors block's step (s)
%
%#codegen

st.sens = sensor_state();
st.gate = gate_state(n_fault_rows, dt);
end

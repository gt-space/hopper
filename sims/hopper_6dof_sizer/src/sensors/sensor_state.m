function st = sensor_state()
% SENSOR_STATE  Initial state for sensor_models: nothing sampled yet.
%
%#codegen

L = sensor_layout();
st.meas  = zeros(L.N_MEAS, 1);
st.truth = zeros(L.N_MEAS, 1);
st.stamp = -ones(L.N_SENSORS, 1);   % < 0: sample on the first frame
end

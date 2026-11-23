function st = sensor_state()
% SENSOR_STATE  Initial state for sensor_models: nothing sampled yet.
%
%#codegen

L = sensor_layout();
st.meas  = zeros(L.N_MEAS, 1);
st.truth = zeros(L.N_MEAS, 1);
st.stamp = -ones(L.N_SENSORS, 1);   % < 0: sample on the first frame
st.w_prev = zeros(3, 1);            % body rates and time of the last frame, for the
st.t_prev = -1;                     % IMU's angular acceleration (< 0: none yet)
end

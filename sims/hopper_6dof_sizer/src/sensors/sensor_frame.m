function [meas, valid, stamp, truth, st] = sensor_frame(x, f_b, dcm_be, cg, noise, t, P, st)
%SENSOR_FRAME  One step of the shared Sensors block: sensors, then faults.
%   [meas, valid, stamp, truth, st] = sensor_frame(x, f_b, dcm_be, cg, noise, t, P, st)
%     x       the controller's 13-state, truth (sim_pos_offset convention)
%     f_b     specific force, body (m/s^2)
%     dcm_be  NED -> body, truth
%     cg      Mass_Properties CG [x; y; z] from the base, body (m)
%     noise   N_MEAS standard-normal samples (seeded Random Number block)
%     P       sensor_params: .S (sim_sensor_params), .faults (fault
%             table), .inject (FAULT_INJECTION on)
%     st      .sens (sensor_state), .gate (gate_state)
%   meas, valid, stamp   what the flight computer sees (sensor_layout),
%                        after the fault table's sensor faults
%   truth                the same channels without noise or faults
%
%   Runs at the fastest sensor rate (the IMU, 500 Hz on the hopper); each
%   sensor samples at its own period and holds between samples. FDIR (100
%   Hz) and the navigation filter (500 Hz) both read this block, so they
%   and the fault injection see the same measurements.
%
%#codegen

x = x(:);
xs.pos_ned = x(1:3); xs.vel_ned = x(4:6); xs.w_b = x(7:9);
xs.dcm_be = dcm_be; xs.f_b = f_b(:); xs.cg = cg(:);
[meas, truth, valid, stamp, st.sens] = sensor_models(xs, noise, t, P.S, st.sens);
if P.inject
    [meas, valid, stamp, st.gate] = sensor_gate(meas, valid, stamp, t, P.faults, st.gate);
end
end

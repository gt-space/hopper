function g = gate_state(n_rows, dt)
% GATE_STATE  Initial state for sensor_gate / actuator_gate.
%
%   n_rows  fault-table rows; dt  gate frame period (s)
%
%#codegen

L = sensor_layout();
g.held   = false(n_rows, 1);
g.value  = zeros(L.N_MEAS, n_rows);
g.stamp  = zeros(L.N_SENSORS, n_rows);
g.held_u = false(n_rows, 1);
g.u_value = zeros(n_rows, 1);
g.u_last  = zeros(n_rows, 1);
g.dt = dt;
end

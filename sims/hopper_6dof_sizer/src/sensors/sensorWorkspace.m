function W = sensorWorkspace(W, variants)
%SENSORWORKSPACE  Add what the shared Sensors block reads to W.
%   W = sensorWorkspace(W, variants)   W from hopperWorkspace, variants from
%                                      variantConfig
%
%   W.sensor_params  .S (sim_sensor_params, with the vehicle's mounting
%                    and IMU rate), .faults (the fault-injection table,
%                    W.faults, default none), .inject (FAULT_INJECTION on),
%                    .dt (the block's step: the IMU period)
%
% The block runs when FDIR is on or NAVIGATION isn't none; FDIR and the
% navigation filter both read it (sensor_frame).

A = variantAxes();
S = sim_sensor_params(W.IN);
if ~isfield(W, 'faults') || isempty(W.faults)
    W.faults = faultScenario();
end
P.S = S;
P.faults = W.faults;
P.inject = A.FAULT_INJECTION(variants.VARIANT_FAULT_INJECTION) == "on";
P.dt = S.imu.period;
if abs(S.frame / P.dt - round(S.frame / P.dt)) > 1e-9 || abs(P.dt / 1e-3 - round(P.dt / 1e-3)) > 1e-9
    error('sensorWorkspace:rate', ['The IMU period (%g s) must divide the FDIR frame (%g s) and be a ' ...
        'multiple of the 1 ms solver step.'], P.dt, S.frame);
end
W.sensor_params = P;
end

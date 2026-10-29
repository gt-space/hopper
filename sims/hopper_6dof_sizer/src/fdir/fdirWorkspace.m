function W = fdirWorkspace(W, variants)
%FDIRWORKSPACE  Add what the FDIR and FAULT_INJECTION choices read to W.
%   W = fdirWorkspace(W, variants)   W from hopperWorkspace, variants from
%                                    variantConfig
%
%   W.fdir_params  sensors, monitor thresholds, mode-manager parameters,
%                  the numeric fault table, the abort profile, and the
%                  fault-injection table (W.faults, default none)
%   W.K_abort      R2/R3 fixed-gain regulator: the nominal LQR gain at
%                  mid-hover (an LQR solution at the operating point the
%                  abort flies)
%   W.K1flat_lo    R1 low-bandwidth gains: half the nominal LQR gains.
%                  LQR keeps its stability margin down to half gain, so
%                  this stays stable; re-solving with a larger R is the
%                  refinement (docs/FDIR_ARCHITECTURE.md).
%   W.fdir_catalog the fault table as read, for reports
%
% Set W.faults = faultScenario(...) before calling buildSimInput to
% inject faults.

A = variantAxes();
IN = W.IN;
S = sim_sensor_params(IN);

P.S = S;
P.mon = gnc_monitor_params(IN, S);
P.modes = fdir_mode_params(IN);
P.modes.contact_alt = S.lidar.mount_x / cosd(S.lidar.cant) + 0.05;   % beam range on the pad + margin
[P.cat, W.fdir_catalog] = fdir_load_catalog();
P.h_cut = P.modes.h_cut;
T = IN.mission.trajectory;
P.abort.v_land = max(T.v_td, 0.1);                         % m/s sink at touchdown
P.abort.v_max  = 0.7 * IN.mission.max_descent_vel;        % m/s fastest abort sink (the nominal descent averages ~4)
P.abort.tau    = 2.5;                                      % s, sink = height / tau
P.abort.a_ramp = 2;                                        % m/s^2: an abort mid-ascent must stop the climb quickly, propellant is short
P.abort.dt     = S.frame;
if ~isfield(W, 'faults') || isempty(W.faults)
    W.faults = faultScenario();
end
P.faults = W.faults;
P.inject = A.FAULT_INJECTION(variants.VARIANT_FAULT_INJECTION) == "on";
P.ref0_err = norm(W.Target_Trajectory(1, 1:3));
P.thrust_max = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(2);
P.thrust_min = IN.propulsion.nominal_thrust * IN.propulsion.throttle_range(1);
P.g = 9.80665;
P.cutoff_alt = IN.mission.engine_cutoff_alt;
W.fdir_params = P;

ih = find(W.tgrid >= T.t_up + T.t_hover / 2, 1);
W.K_abort = reshape(W.K1flat(ih, :), 13, 4)';
W.K1flat_lo = 0.5 * W.K1flat;
end

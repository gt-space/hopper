function REQ = tvcRequirements()
%TVCREQUIREMENTS  GNC-derived gimbal requirements handed to the TVC team.
%   Each value carries its derivation; tvcTrade checks both actuators
%   against these. Update here, not in the trade script.

% ---- Loop / vehicle numbers (GNC PDR) ----
REQ.wc_rad_s  = 12.25;     % pitch-loop crossover (PDR slide 164, 1.95 Hz)
REQ.wn_act    = 36;        % rad/s, actuator small-signal bandwidth the LQT was tuned with
REQ.axial_g   = 1.6;       % peak axial accel (g) for the engine-weight hinge term (TBR)
l_gimbal_cg   = 1.040 - 0.45;   % m, dry CG to gimbal (hopper_config; worst = shortest arm)

% ---- Range ----
% trim for the 2.5 cm stress thrust-line offset + 5 mm lateral CG, plus the
% 6DOF/MC control demand (<~1.5 deg, PDR slide 164) and wind/slosh margin.
REQ.trim_deg  = rad2deg(atan((0.025 + 0.005) / l_gimbal_cg));   % ~2.9 deg
REQ.ctrl_deg  = 3.0;       % control authority above trim (2x observed MC peak)
REQ.range_deg = 8;         % commanded +/- per axis, both axes simultaneously (square)
REQ.hardstop_deg = 9.5;    % mechanical stop: above range + overshoot, bounds a hard-over

% ---- Rate ----
% stay linear (no rate saturation) up to A_lin = 2 deg at crossover:
% rate >= A_lin * wc = 2 deg * 12.25 rad/s = 24.5 deg/s  -> 25 deg/s at load,
% each axis, both axes moving together.
REQ.A_lin_deg = 2;
REQ.rate_dps  = 25;

% ---- Acceleration ----
% follow the same linear amplitude as the rate requirement at crossover:
% A_lin * wc^2 = 2 deg * 12.25^2 = 300 deg/s^2.  Was 900 (25 deg/s * 36
% rad/s, full rate within one actuator time constant), which the TVC team
% cannot reach and the vehicle does not need: with the drive's acceleration
% ramp at 150-900 deg/s^2 (AG 9113 recommended, 120 in-loop runs, results/
% tvc_accel_loop.mat) landing is unchanged to 1 mm in every MC and stress
% case; only transient tracking error grows (2.5 cm offset worst 0.43 deg
% unlimited, 0.87 at 300). Below ~200 the 2 Hz response to a 2 deg swing
% collapses (phase past -150 deg), so 300 is the floor, not 150.
% Set the drive's acceleration ramp no higher than the demonstrated value.
REQ.accel_dps2 = 300;

% ---- Response / resolution / load ----
REQ.rise05_ms    = 65;     % 0.5 deg step, 10-90 %, latency excluded (63 ms at the 300 deg/s^2 ramp; was 60)
REQ.phase2Hz_deg = -40;    % actuator phase at 2 Hz, sine amplitude <= 1 deg, no worse than this (current 6DOF -56)
                           % (flight swings at crossover are well under 1 deg: peak in-loop rate 6.7 deg/s)
REQ.res_deg      = 0.05;   % <= half the 6DOF 0.1 deg quantizer; 4x below deadband 0.2
REQ.force_margin = 1.5;    % rated force / worst hinge load (2.5 cm stress offset, 8 deg, 1.1x thrust)
REQ.force_margin_nom = 2.5;% rated force / nominal hinge load (1 cm offset)
REQ.latency_ms   = 10;     % command-to-motion, Ethernet + servo (TBR with E+L)
end

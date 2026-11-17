function P = tvcActuatorParams(name, opts)
%TVCACTUATORPARAMS  Erhardt+Leimer AG 91.3 / 93.3 TVC actuator + gimbal values.
%   P = tvcActuatorParams('AG9113')   PDR baseline (1000 N, 60 mm/s)
%   P = tvcActuatorParams('AG9313')   sizer's current pick (3000 N, 30 mm/s)
%   P = tvcActuatorParams(..., RodEndZ=0.07)  move the rod end (arm study)
%
%   Sources
%     Datasheet DAT--251370-EN-10 / manual BEA--251942-EN-05: Vbus, Imax,
%       Fmax, rated speed, stroke, mid-length M, accuracy (<0.05 mm), mass
%     GNC PDR slides 165-172: model form, gimbal/hinge loads, servo, latency
%     HopperTVCActuatorSizer.m: tail-stock / rod-end coordinates (in)
%   Derived (same method as PDR slide 166 for the AG 9113):
%     no-load speed  v0 = (7/6) * rated speed   (PDR estimate 70 vs 60 mm/s)
%     Ke = Vbus / v0,  KF = Fmax / Imax,  R0 = (Vbus - Ke*v_rated) / Imax
%   Estimated (TBR by bench test): L, rod-referred rotor mass m_r,
%     Coulomb friction, backlash, servo gains (scaled for same rise time).

arguments
    name (1,:) char {mustBeMember(name, {'AG9113', 'AG9313'})}
    opts.RodEndZ double = []      % m, rod-end radial offset from engine axis
    opts.TailStockX double = []   % m, tail-stock axial station
    opts.FitLength (1,1) logical = false   % set L0 = datasheet mid length M
    opts.Layout (1,:) char {mustBeMember(opts.Layout, {'sizer', 'recommended'})} = 'sizer'
end
in = 0.0254;
P.name = name;
P.layout = opts.Layout;
if strcmp(opts.Layout, 'recommended')
    % tvcTrade 2026-11-15: arm that balances rate and force margin, with
    % the actuator at mid-stroke at 0 deg (r_eff 122.5 / 57.5 mm)
    REC = struct('AG9113', 0.1113, 'AG9313', 0.0302);   % m, rod-end radius
    if isempty(opts.RodEndZ), opts.RodEndZ = REC.(name); end
    opts.FitLength = true;
end

% ---- Datasheet ----
P.Vbus = 24;                        % V
P.q_mech = 0.05e-3;                 % m, mechanical position accuracy
switch name
    case 'AG9113'
        P.Imax = 5.6;  P.Fmax = 1000;  P.v_rated = 0.060;
        P.stroke = 0.050;           % m, +/- nominal actuating travel
        P.M = 0.360;                % m, mid-position pin-to-pin length
        P.mass = 3.9;               % kg (+/-50 mm)
        P.m_r = 25;                 % kg, rod-referred rotor+screw (PDR est.)
        % Sizer geometry (commented-out 9113 block), in
        RE = [-3.35; 0; 5.00];  TS = [8.80; 0; 6.45];
    case 'AG9313'
        P.Imax = 7.7;  P.Fmax = 3000;  P.v_rated = 0.030;
        P.stroke = 0.050;
        P.M = 0.390;
        P.mass = 5.0;
        % Half the speed for the same motor class -> ~half the lead, and the
        % reflected inertia goes as 1/lead^2 -> ~4x the AG 9113 value.
        P.m_r = 100;                % kg, ESTIMATED, sensitivity 50-200 kg
        RE = [-3.35; 0; 3.50];  TS = [11.75; 0; 6.45];
end
P.v0 = 7/6 * P.v_rated;                       % m/s no-load (est.)
P.Ke = P.Vbus / P.v0;                         % V s/m
P.KF = P.Fmax / P.Imax;                       % N/A
P.R0 = (P.Vbus - P.Ke * P.v_rated) / P.Imax;  % ohm
P.eta = P.KF / P.Ke;
P.L = 1e-3;                                   % H (est.)
P.Fc = 20 * P.Fmax / 1000;                    % N Coulomb (est., scales w/ size)
P.backlash = deg2rad(0.03);                   % rad (est., rod ends + pins)

% ---- Geometry (engine frame, origin at gimbal pivot, x along engine) ----
P.RE_P = RE * in;  P.TS_P = TS * in;          % pitch actuator (x-z plane)
if ~isempty(opts.RodEndZ),    P.RE_P(3) = opts.RodEndZ; end
if ~isempty(opts.TailStockX), P.TS_P(1) = opts.TailStockX; end
if opts.FitLength
    % move the tail stock axially so the neutral pin-to-pin length is the
    % datasheet mid-position length M (actuator at mid-stroke at 0 deg)
    P.TS_P(1) = P.RE_P(1) + sqrt(P.M^2 - (P.TS_P(3) - P.RE_P(3))^2);
end
P.RE_Y = [P.RE_P(1); P.RE_P(3); 0];           % yaw actuator (x-y plane)
P.TS_Y = [P.TS_P(1); P.TS_P(3); 0];
P.L0 = norm(P.RE_P - P.TS_P);                 % m, neutral pin-to-pin length
u = (P.RE_P - P.TS_P) / P.L0;
P.r_eff = abs(P.RE_P(3) * u(1) - P.RE_P(1) * u(3));   % m, dL/d(delta) at neutral

% ---- Gimbal + engine loads (PDR slides 166, 172) ----
P.Jg   = 0.042;            % kg m^2, engine about pivot (CAD tent.)
P.kh   = 49.8;             % N m/rad, LOx feedline
P.tauf = 2.82;             % N m, gimbal pin friction
P.me = 2.35;  P.rcg = 0.090;   % kg, m
P.thrust = 1.1 * 500 * 4.44822;   % N, max throttle
P.e_nom = 0.010;  P.e_stress = 0.025;   % m, thrust-line offset

% ---- Built-in servo (1 kHz) and command path ----
P.dt_c = 1e-3;             % s, controller period
P.lat  = 0.010;            % s, command latency (sim value, TBR with E+L)
% PDR gains for the AG 9113 at r_eff = 0.114 m, m_r = 25 kg; scale with the
% plant gain (m_r r_eff / KF) so the small-signal loop shape is kept.
g = (P.m_r * P.r_eff / P.KF) / (25 * 0.114 / 178);
P.Kp = 4000 * g;  P.Ki = 20000 * g;  P.Kd = 60 * g;    % V/rad, V/(rad s), V s/rad

% ---- Capability summary (closed form) ----
% slew at load torque T: i = T/(r KF); v = (Vbus - R0 i)/Ke; rate = v/r
P.rate_noload = P.v0 / P.r_eff;                               % rad/s
P.torque_max  = P.Fmax * P.r_eff;                              % N m (rated)
P.res_angle   = P.q_mech / P.r_eff;                            % rad
P.range_stroke = P.stroke / P.r_eff;                           % rad (linearized)
end

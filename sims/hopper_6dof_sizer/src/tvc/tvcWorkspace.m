function W = tvcWorkspace(W, variants, opts)
%TVCWORKSPACE  TVCEM for the TVC_MODEL variant (Actuators_TVC_RCS/TVC_Actuator).
%   W = tvcWorkspace(W, variants) adds W.TVCEM, the electromechanical
%   actuator parameters, unless W already has one (a caller such as
%   tvcTrade can set its own actuator or layout first).
%     TVC_MODEL = second_order   the PDR 6DOF chain; TVCEM is a placeholder
%                 ag9113, ag9313 Erhardt+Leimer AG 91.3 / 93.3 on the
%                                recommended layout (tvcActuatorParams)
%   Latency comes from IN.tvc.actuator_latency (dispersed in Monte Carlo);
%   the actuator load uses the thrust-line offset IN.off_axis_* with a
%   1 cm floor (the PDR actuator load case), and the gimbal range is
%   the GNC requirement, +/-8 deg.

%   W = tvcWorkspace(W, variants, Layout='sizer', Offset=0.025) replaces
%   any TVCEM in W: the actuator's layout and a fixed actuator load offset
%   (m) for trade studies.  AccelLimit=deg2rad(300) sets the drive's
%   acceleration ramp (rad/s^2; default Inf, no ramp).
arguments
    W struct
    variants struct
    opts.Layout (1,:) char = ''
    opts.Offset double = []
    opts.AccelLimit (1,1) double = Inf
end
if isfield(W, 'TVCEM') && isempty(opts.Layout) && isempty(opts.Offset) && isinf(opts.AccelLimit), return; end
A = variantAxes();
choice = "second_order";
if isfield(variants, 'VARIANT_TVC_MODEL')
    choice = A.TVC_MODEL(variants.VARIANT_TVC_MODEL);
end
name = 'AG9113';
if choice == "ag9313", name = 'AG9313'; end
layout = opts.Layout;  if isempty(layout), layout = 'recommended'; end
P = tvcActuatorParams(name, Layout=layout);

lat = P.lat;  e = [0.01; 0.01];
if isfield(W, 'IN')
    if isfield(W.IN, 'tvc') && isfield(W.IN.tvc, 'actuator_latency')
        lat = W.IN.tvc.actuator_latency;
    end
    if isfield(W.IN, 'off_axis_x'), e(1) = max(abs(W.IN.off_axis_x), 0.01); end   % body z offset -> pitch
    if isfield(W.IN, 'off_axis_y'), e(2) = max(abs(W.IN.off_axis_y), 0.01); end   % body y offset -> yaw
end
if ~isempty(opts.Offset), e = [opts.Offset; opts.Offset]; end
W.TVCEM = tvcEmParams(P, Latency=lat, Offset=e, AccelLimit=opts.AccelLimit);
% the Electromech choice delays its inputs one servo period (to break the
% algebraic loop through thrust); take it out of the delay line
W.TVCEM.lat_steps = max(W.TVCEM.lat_steps - 1, 0);
end

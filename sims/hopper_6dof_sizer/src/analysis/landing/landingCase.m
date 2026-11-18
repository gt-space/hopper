function C = landingCase(opts)
%LANDINGCASE  A touchdown state for landingSim.
%   C = landingCase(SinkRate=1.5, Lateral=0.5, Heading='leg', Tilt=2)
%
%   Heading and TiltToward are body azimuths about x from +y (deg), or
%   'leg' / 'between' for the azimuth of leg 1 or halfway to leg 2:
%     'leg'      one foot leads (contact order 1-2-1 in the Ansys naming)
%     'between'  two feet lead (2-2)
%   Attitude sets which side is low at touchdown (that foot lands first):
%     'braking'  top tilted away from the drift, as when the vehicle is
%                killing its sideways speed: the trailing foot lands first
%                and the vehicle rocks forward onto the leading foot
%     'leading'  top tilted into the drift (e.g. after a gust): the leading
%                foot lands first and digs in, the tip-over direction
%   TiltToward (azimuth of the low side) overrides Attitude. Mass properties
%   default to the 145 kg design mass of the leg hand calcs, with the dry CG
%   and inertia scaled from the CAD.

arguments
    opts.Name (1,:) char = ''
    opts.SinkRate (1,1) double = 1.5          % m/s, vertical, positive down
    opts.Lateral (1,1) double = 0             % m/s, horizontal
    opts.Heading = 'leg'
    opts.Tilt (1,1) double = 0                % deg, body axis from vertical
    opts.Attitude (1,:) char {mustBeMember(opts.Attitude, {'braking', 'leading'})} = 'braking'
    opts.TiltToward = []
    opts.Rate (3,1) double = [0; 0; 0]        % rad/s, body
    opts.Mass (1,1) double = 145              % kg
    opts.CgX double = []                      % m above the pad bottoms
    opts.Inertia double = []                  % 3x3 kg m^2, body
    opts.Thrust (1,1) double = 0              % N along body x at contact
    opts.ThrustTau (1,1) double = 0.05        % s, thrust decay after contact
    opts.ThrustProfile = []                   % [t; fraction]: thrust after the cutoff command, as a fraction
                                              % of Thrust (then Thrust is the thrust at the command); replaces ThrustTau
    opts.ThrustTime0 (1,1) double = 0         % s since the cutoff command at contact (< 0: cut after contact)
    opts.Vehicle = []
end
VC = opts.Vehicle;
if isempty(VC), VC = vehicleConfig('hopper'); end
leg_az = atan2d(VC.legs.geom.node(3), VC.legs.geom.node(2));

C.sink = opts.SinkRate;
C.lateral = opts.Lateral;
C.heading = azimuth(opts.Heading, leg_az);
C.attitude = opts.Attitude;
if ~isempty(opts.TiltToward)
    C.tilt_toward = azimuth(opts.TiltToward, leg_az);
    C.attitude = 'custom';
elseif strcmp(opts.Attitude, 'braking')
    C.tilt_toward = C.heading + 180;            % trailing side low
else
    C.tilt_toward = C.heading;                  % leading side low
end
C.tilt = opts.Tilt;
C.rate = opts.Rate;
C.mass = opts.Mass;
C.cg_x = opts.CgX;
if isempty(C.cg_x), C.cg_x = VC.mass.dry_cg; end
C.I = opts.Inertia;
if isempty(C.I), C.I = VC.mass.dry_I * opts.Mass / VC.mass.dry; end
C.thrust = opts.Thrust;
C.thrust_tau = opts.ThrustTau;
C.thrust_profile = opts.ThrustProfile;
C.thrust_t0 = opts.ThrustTime0;
C.name = opts.Name;
if isempty(C.name)
    C.name = sprintf('v%.2f_h%.2f_%s_t%g_%s', C.sink, C.lateral, headingName(opts.Heading), C.tilt, C.attitude);
end
end

function a = azimuth(h, leg_az)
if ischar(h) || isstring(h)
    switch char(h)
        case 'leg',     a = leg_az;
        case 'between', a = leg_az + 45;
        otherwise, error('landingCase:heading', 'Heading must be ''leg'', ''between'' or degrees.');
    end
else
    a = h;
end
end

function s = headingName(h)
if ischar(h) || isstring(h), s = char(h); else, s = sprintf('%gdeg', h); end
end

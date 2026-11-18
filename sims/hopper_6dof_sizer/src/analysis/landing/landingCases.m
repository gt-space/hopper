function cases = landingCases(opts)
%LANDINGCASES  Touchdown case matrix for the landing-loads sweep.
%   cases = landingCases()   struct array with fields case (landingCase)
%                            and set ('nominal' | 'offnominal')
%
%   nominal:    the GNC touchdown requirement (mission_inputs.m): sink
%               <= 1.5 m/s, lateral <= 0.5 m/s, tilt <= 2 deg
%   offnominal: the leg hand calcs' crush-core case, 10 ft/s sink and
%               2 ft/s lateral, with the same attitudes and a 5 deg tilt
%
%   Heading 'leg': one foot leads (Ansys 1-2-1); 'between': two lead (2-2).
%   Each drift case runs in both attitudes (landingCase): 'braking', top
%   tilted away from the drift so the trailing foot lands first, and
%   'leading', top tilted into the drift so the leading foot lands first.

arguments
    opts.Mass (1,1) double = 145
    opts.Vehicle = []
end
VC = opts.Vehicle;
if isempty(VC), VC = vehicleConfig('hopper'); end
m = opts.Mass;
mk = @(name, set, varargin) struct('case', landingCase('Name', name, 'Mass', m, 'Vehicle', VC, varargin{:}), 'set', set);

cases = [ ...
    mk('N1 flat',                          'nominal',    SinkRate=1.5)
    mk('N2 drift to leg, braking',         'nominal',    SinkRate=1.5, Lateral=0.5, Heading='leg', Tilt=2, Attitude='braking')
    mk('N2 drift to leg, leading',         'nominal',    SinkRate=1.5, Lateral=0.5, Heading='leg', Tilt=2, Attitude='leading')
    mk('N3 drift between, braking',        'nominal',    SinkRate=1.5, Lateral=0.5, Heading='between', Tilt=2, Attitude='braking')
    mk('N3 drift between, leading',        'nominal',    SinkRate=1.5, Lateral=0.5, Heading='between', Tilt=2, Attitude='leading')
    mk('O1 flat',                          'offnominal', SinkRate=3.048)
    mk('O2 drift to leg, braking',         'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='leg', Tilt=2, Attitude='braking')
    mk('O2 drift to leg, leading',         'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='leg', Tilt=2, Attitude='leading')
    mk('O3 drift between, braking',        'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='between', Tilt=2, Attitude='braking')
    mk('O3 drift between, leading',        'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='between', Tilt=2, Attitude='leading')
    mk('O4 drift to leg 5deg, braking',    'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='leg', Tilt=5, Attitude='braking')
    mk('O4 drift to leg 5deg, leading',    'offnominal', SinkRate=3.048, Lateral=0.6096, Heading='leg', Tilt=5, Attitude='leading')];
end

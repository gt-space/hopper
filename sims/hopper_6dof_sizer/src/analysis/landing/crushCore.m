function C = crushCore(placement, opts)
%CRUSHCORE  A crush-core configuration for landingGear.
%   C = crushCore('none')
%   C = crushCore('main', Plateau=15e3)            inline in each main strut
%   C = crushCore('pad',  Plateau=8e3, Stroke=0.05) under each footpad
%
%   Force law (compression, per core): the core holds elastically until
%   the load reaches Peak, falls linearly to Plateau over the first Trigger
%   of crush, stays at Plateau to the usable stroke (Densify x Stroke), then
%   stiffens like a solid. Crush is permanent; unloading is elastic.
%
%   'main': in series with the main strut. Aux struts run in tension on
%           landing, so a compression core there would not crush.
%   'pad':  between the foot fitting and the ground, normal to the ground.
%
%   Stroke defaults to Inf so a sweep reports the stroke a load level needs
%   (stroke used / Densify) instead of bottoming out.
%
%   Before it crushes the core is a spring, Modulus * A / ElasticLength,
%   with A = Plateau / CrushStress (the core area that gives this plateau).
%   Defaults are aluminium honeycomb (5052, ~3 pcf: 1.8 MPa crush, ~0.5 GPa
%   stabilized compressive modulus) and a 0.13 m core. The spring is in
%   series with the strut (main) or the ground contact (pad), so it softens
%   the foot slap. Stiffness overrides it; Inf makes the core rigid.

arguments
    placement (1,:) char {mustBeMember(placement, {'none', 'main', 'pad'})}
    opts.Plateau (1,1) double {mustBeNonnegative} = 0      % N
    opts.PeakRatio (1,1) double {mustBeGreaterThanOrEqual(opts.PeakRatio, 1)} = 1.0   % peak / plateau (1 = pre-crushed)
    opts.Trigger (1,1) double {mustBePositive} = 2e-3      % m, crush over which the peak decays
    opts.Stroke (1,1) double {mustBePositive} = Inf        % m, core length
    opts.Densify (1,1) double {mustBePositive} = 0.75      % usable fraction of the length
    opts.DensifyStiffness (1,1) double {mustBePositive} = 5e7   % N/m past densification
    opts.CrushStress (1,1) double {mustBePositive} = 1.8e6      % Pa, core plateau stress (sets the core area)
    opts.Modulus (1,1) double {mustBePositive} = 517e6          % Pa, core elastic modulus
    opts.ElasticLength (1,1) double {mustBePositive} = 0.13     % m, core length for the elastic spring
    opts.Stiffness double = []                                  % N/m, overrides the computed spring
end
C.placement = placement;
C.plateau = opts.Plateau;
C.peak = opts.Plateau * opts.PeakRatio;
C.trigger = opts.Trigger;
C.length = opts.Stroke;
C.s_max = opts.Densify * opts.Stroke;
C.densify = opts.Densify;
C.k_dens = opts.DensifyStiffness;
C.area = opts.Plateau / opts.CrushStress;
if ~isempty(opts.Stiffness)
    C.k_elastic = opts.Stiffness;
elseif strcmp(placement, 'none')
    C.k_elastic = Inf;
else
    C.k_elastic = opts.Modulus * C.area / opts.ElasticLength;
end
if ~strcmp(placement, 'none') && C.plateau <= 0
    error('crushCore:plateau', 'A %s core needs a positive Plateau force.', placement);
end
end

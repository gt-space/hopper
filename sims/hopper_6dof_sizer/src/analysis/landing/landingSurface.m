function S = landingSurface(name, opts)
%LANDINGSURFACE  Ground model presets for landingSim.
%   S = landingSurface('concrete')   pad on concrete
%   S = landingSurface('soil_firm')  compacted soil / dense sand
%   S = landingSurface('soil_soft')  loose sand / soft soil
%   S = landingSurface('rigid')      near-rigid ground, to match an FEA drop
%   Any field can be overridden: landingSurface('concrete', mu=0.4)
%
%   'linear' surfaces: normal force k*pen + c*pen_rate, never pulling.
%   'soil' surfaces: elastic-plastic. Loading follows the backbone
%   ks*A_pad*pen (modulus of subgrade reaction times pad area); unloading
%   is Ku_ratio times stiffer, and the pad leaves a permanent dent.
%   Friction: a stick-slip spring (kt) capped at mu * normal force; soil's
%   mu stands in for the pad ploughing into the ground.
%
%   Values are textbook ranges, not site measurements:
%     concrete k: pad and fitting local compliance dominate (2e8 N/m)
%     subgrade ks: loose sand 5-16, medium 10-80, dense 64-128 MN/m^3

arguments
    name (1,:) char {mustBeMember(name, {'rigid', 'concrete', 'soil_firm', 'soil_soft'})}
    opts.k double = []
    opts.zeta double = []
    opts.mu double = []
    opts.kt double = []
    opts.ks double = []
    opts.Ku_ratio double = []
end
A_pad = pi/4 * 0.175^2;
switch name
    case 'rigid'
        S = struct('type', 'linear', 'k', 1e9, 'zeta', 0.05, 'mu', 0.3, 'kt', 1e8, 'ks', 0, 'Ku_ratio', 1);
    case 'concrete'
        S = struct('type', 'linear', 'k', 2e8, 'zeta', 0.2, 'mu', 0.6, 'kt', 5e7, 'ks', 0, 'Ku_ratio', 1);
    case 'soil_firm'
        S = struct('type', 'soil', 'k', 0, 'zeta', 0.2, 'mu', 0.8, 'kt', 5e6, 'ks', 100e6, 'Ku_ratio', 10);
    case 'soil_soft'
        S = struct('type', 'soil', 'k', 0, 'zeta', 0.2, 'mu', 0.8, 'kt', 1e6, 'ks', 15e6, 'Ku_ratio', 10);
end
S.name = name;
S.A_pad = A_pad;
for f = ["k", "zeta", "mu", "kt", "ks", "Ku_ratio"]
    if ~isempty(opts.(f)), S.(f) = opts.(f); end
end
if strcmp(S.type, 'soil')
    S.kbb = S.ks * A_pad;                % N/m, loading backbone per pad
    S.Ku = S.Ku_ratio * S.kbb;           % N/m, unloading
else
    S.kbb = S.k;  S.Ku = S.k;
end
end

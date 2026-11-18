function G = landingGear(opts)
%LANDINGGEAR  The as-drawn landing gear for the landing-loads model.
%   G = landingGear()                     pinned struts, no crush cores
%   G = landingGear(Core=crushCore('main', Plateau=15e3))
%   G = landingGear(Joints='fixed')       struts clamped at both ends (frame)
%
%   Each leg is a foot node (pad, fitting and half the strut mass) held to
%   the body by three struts (main, auxA, auxB) from VC.legs.geom. Strut
%   lengths are exact (large displacement), so a stroking crush core moves
%   the foot along the real mechanism.
%
%   Joints:
%     'pinned'  struts carry axial load only (clevis pins at both ends, as
%               in the hand calcs). The three struts lie within 10 deg of
%               each other, so a pinned leg is soft across the strut bundle
%               (about 0.9 MN/m vertical) and the main strut sees about
%               5x the vertical foot load.
%     'fixed'   struts also bend, clamped at the body and at a foot fitting
%               that is free to rotate (ball-jointed pad). A stiff bound.
%
%   G fields (body axes, m, from the footpad bottoms):
%     node, pad, upper (3x3xn), L0, u0, k, c, A, I, Pcr, m_foot, Kb, Mmap

arguments
    opts.Vehicle = []
    opts.Joints (1,:) char {mustBeMember(opts.Joints, {'pinned', 'fixed'})} = 'pinned'
    opts.Core struct = crushCore('none')
    opts.StrutDamping (1,1) double {mustBeNonnegative} = 0.02   % damping ratio of each strut on the foot mass
    opts.FittingMass double = []                                % kg per foot; default VC.legs.foot_fitting_mass
    opts.FootDamping (1,1) double {mustBeNonnegative} = 0.05    % damping ratio of the foot on the leg's soft mode (joint friction)
end

VC = opts.Vehicle;
if isempty(VC), VC = vehicleConfig('hopper'); end
geo = VC.legs.geom;
st = VC.legs.strut;
n = numel(geo.turns);

G.n = n;
G.joints = opts.Joints;
G.core = opts.Core;
G.E = st.E;  G.sy = st.sy;
G.A = pi/4 * (st.od.^2 - st.id.^2);           % 1x3
G.I = pi/64 * (st.od.^4 - st.id.^4);          % 1x3, area moment (not polar)
G.c_out = st.od / 2;
G.pad_diameter = geo.pad_diameter;
G.names = {'main', 'auxA', 'auxB'};

G.node = zeros(3, n);  G.pad = zeros(3, n);  G.upper = zeros(3, 3, n);
for k = 1:n
    Rk = rotx_(geo.turns(k));
    G.node(:, k) = Rk * geo.node;
    G.pad(:, k)  = Rk * geo.pad;
    G.upper(:, :, k) = Rk * geo.upper;
end
d = reshape(G.node, 3, 1, n) - G.upper;       % column i = node - upper_i
G.L0 = squeeze(vecnorm(d, 2, 1));             % 3xn
G.u0 = d ./ reshape(G.L0, 1, 3, n);           % unit, upper -> node
G.k  = (G.E * G.A.') ./ G.L0;                 % 3xn, N/m
G.Pcr = pi^2 * G.E * G.I.' ./ G.L0.^2;        % 3xn, pinned-pinned Euler load

strut_mass = st.rho * G.A * G.L0(:, 1);       % kg per leg
fit = opts.FittingMass;
if isempty(fit), fit = VC.legs.foot_fitting_mass; end
G.strut_mass = strut_mass;
G.m_foot = fit + strut_mass / 2;
G.c = 2 * opts.StrutDamping * sqrt(G.k * G.m_foot);
% The pinned leg is nearly a mechanism across the strut bundle; strut
% damping does not reach that mode, so the foot also gets an isotropic
% damper relative to the body (joint friction, ASSUMPTION).
u1 = G.u0(:, :, 1);
K1 = u1 * diag(G.k(:, 1)) * u1.';
G.c_foot = 2 * opts.FootDamping * sqrt(min(eig(K1)) * G.m_foot);

% Bending (fixed joints): condensed 3x3 translational stiffness of the
% foot node without the axial terms, and maps from foot displacement to
% the bending moment at each strut end.
G.Kb = zeros(3, 3, n);
G.Mmap = zeros(2, 3, 2, 3, n);                % (moment comp, disp, end, strut, leg)
if strcmp(opts.Joints, 'fixed')
    Gs = G.E / 2.6;  J = 2 * G.I;
    for k = 1:n
        K66 = zeros(6);  Ke = cell(1, 3);  T = cell(1, 3);
        for i = 1:3
            [kl, lam] = beamNoAxial(G.E, G.I(i), Gs, J(i), G.L0(i, k), -G.u0(:, i, k));
            Tb = blkdiag(lam, lam, lam, lam);
            Ke{i} = kl;  T{i} = Tb;
            Kg = Tb.' * kl * Tb;
            K66 = K66 + Kg(1:6, 1:6);         % node is end 1; the upper end is clamped
        end
        Ktt = K66(1:3, 1:3); Ktr = K66(1:3, 4:6); Krr = K66(4:6, 4:6);
        G.Kb(:, :, k) = Ktt - Ktr / Krr * Ktr.';
        Rot = -Krr \ Ktr.';                   % node rotation per unit node displacement
        for i = 1:3
            D = [eye(3); Rot; zeros(6, 3)];   % global end DOFs per unit displacement
            f = Ke{i} * T{i} * D;             % local end forces
            G.Mmap(:, :, 1, i, k) = f([5 6], :);
            G.Mmap(:, :, 2, i, k) = f([11 12], :);
        end
    end
end

% Static tip-over angle about the edge between adjacent pads, dry CG
pads = G.pad(2:3, :);
edge_r = min(vecnorm((pads + circshift(pads, -1, 2)) / 2));
G.cg_x = VC.mass.dry_cg;
G.tip_angle_deg = atand(edge_r / G.cg_x);
end

function R = rotx_(deg)
c = cosd(deg); s = sind(deg);
R = [1 0 0; 0 c -s; 0 s c];
end

function [k, lam] = beamNoAxial(E, I, G, J, L, ex)
% 12x12 Euler-Bernoulli beam stiffness, local axes, axial term removed.
% End 1 at the foot node, local x from end 1 to end 2 (= ex).
ref = [0; 0; 1];
if abs(dot(ref, ex)) > 0.9, ref = [0; 1; 0]; end
ey = cross(ref, ex); ey = ey / norm(ey);
ez = cross(ex, ey);
lam = [ex.'; ey.'; ez.'];
k = zeros(12);
a = 12*E*I/L^3; b = 6*E*I/L^2; c4 = 4*E*I/L; c2 = 2*E*I/L; t = G*J/L;
% bending in local x-y (v, rz): DOFs 2,6,8,12
idx = [2 6 8 12];
kb = [a b -a b; b c4 -b c2; -a -b a -b; b c2 -b c4];
k(idx, idx) = kb;
% bending in local x-z (w, ry): DOFs 3,5,9,11 (sign flips on the coupling)
idx = [3 5 9 11];
kb = [a -b -a -b; -b c4 b c2; -a b a b; -b c2 b c4];
k(idx, idx) = kb;
% torsion: DOFs 4,10
k([4 10], [4 10]) = t * [1 -1; -1 1];
end

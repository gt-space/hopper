function R = landingSim(G, S, C, opts)
%LANDINGSIM  Touchdown of the hopper on its as-drawn legs: loads and stroke.
%   R = landingSim(landingGear(), landingSurface('concrete'), landingCase(SinkRate=1.5))
%
%   A rigid body (the vehicle less its feet) on four legs. Each foot node
%   is a point mass joined to the body by the main and two aux struts
%   (exact lengths, axial spring + damping, optional crush core), and to
%   the ground by the surface model (normal spring or elastic-plastic soil,
%   optional crush pad, stick-slip friction). Fixed-step RK4; the plastic
%   states (core strokes, soil dents, friction anchors) are updated after
%   each step by return mapping.
%
%   Inertial frame: X, Y level, Z up, ground at Z = 0. Body x is the thrust
%   axis. Ground force acts at the foot node (pad pivot assumed).
%
%   R.peak  per-strut peak compression / tension (N, compression +), peak
%           bending moment and combined stress (fixed joints), per-foot peak
%           normal and friction force, total ground force, crush stroke,
%           soil dent, peak CG acceleration (g), max tilt, energy absorbed
%   R.hist  decimated time history (opts.Log)
%
%   Raw peaks include the foot-slap and node-rattle spikes (kHz); acc_lp_g,
%   Fn_lp and Fn_total_lp are the same signals through a first-order
%   low-pass at FilterHz, the values to quote for equipment and for a
%   quasi-static ground load.
%
%   The run goes to tEnd (2.5 s, long enough to see a slow tip-over) or
%   stops once the body has been at rest for SettleTime, or has tipped
%   past 80 deg.
%
%   opts: dt (2e-5 s), tEnd (2.5 s), SettleTime (0.1 s), LogEvery (10
%         steps), Log (true), FilterHz (100)

arguments
    G struct
    S struct
    C struct
    opts.dt (1,1) double {mustBePositive} = 2e-5
    opts.tEnd (1,1) double {mustBePositive} = 2.5
    opts.SettleTime (1,1) double {mustBePositive} = 0.1   % s at rest before stopping early
    opts.LogEvery (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.Log (1,1) logical = true
    opts.FilterHz (1,1) double {mustBePositive} = 100   % low-pass for the reported CG g and ground force
end

P = setup(G, S, C);
dt = opts.dt;
nStep = ceil(opts.tEnd / dt);

% ---- Initial state ----
R0 = [0 1 0; 0 0 1; 1 0 0];                  % upright: body x->Z, y->X, z->Y
d_b = [0; cosd(C.tilt_toward); sind(C.tilt_toward)];
ax = cross([1; 0; 0], d_b);  ax = ax / norm(ax);   % top tips toward d_b, so the d_b-side foot is low
Rbi = R0 * axang(ax, deg2rad(C.tilt));
q = dcm2q(Rbi);
padB = P.pad_b;                               % pads relative to CG, body
padI = Rbi * padB;
p = [0; 0; -min(padI(3, :)) + 1e-4];
v = [C.lateral * cosd(C.heading); C.lateral * sind(C.heading); -C.sink];
w = C.rate;
f = p + Rbi * P.node_b;
fv = v + Rbi * cross(repmat(w, 1, 4), P.node_b);
x = [p; v; q; w; f(:); fv(:)];

Z.s  = zeros(3, 4);                           % strut core crush (m)
Z.sp = zeros(1, 4);                           % pad core crush (m)
Z.dp = zeros(1, 4);                           % soil plastic dent (m)
Z.a  = f(1:2, :) + Rbi(1:2, :) * P.off_b;     % friction anchors = contact points

% ---- Peaks and log ----
pk.Nc = zeros(3, 4);  pk.Nt = zeros(3, 4);  pk.M = zeros(3, 4);  pk.sig = zeros(3, 4);
pk.Fn = zeros(1, 4);  pk.Ft = zeros(1, 4);  pk.Fn_total = 0;  pk.acc_g = 0;  pk.acc_axial_g = 0;
pk.tilt_deg = C.tilt;  pk.E_core = 0;  pk.E_pad = 0;
pk.acc_lp_g = 0;  pk.Fn_lp = zeros(1, 4);  pk.Fn_total_lp = 0;
a_lp = zeros(3, 1);  fn_lp = zeros(1, 4);
alpha = dt / (dt + 1 / (2 * pi * opts.FilterHz));   % first-order low-pass
nLog = floor(nStep / opts.LogEvery) + 1;
if opts.Log
    H.t = zeros(1, nLog);  H.z = H.t;  H.vz = H.t;  H.tilt = H.t;  H.acc_g = H.t;  H.acc_lp_g = H.t;
    H.N = zeros(3, 4, nLog);  H.Fn = zeros(4, nLog);  H.Ft = zeros(4, nLog);
    H.s = zeros(3, 4, nLog);  H.sp = zeros(4, nLog);  H.dp = zeros(4, nLog);
end
iLog = 0;
E0 = energy(x, P);

t = 0;  settled = 0;
for n = 0:nStep
    [k1, I1] = rates(x, Z, t, P);
    % peaks from the start-of-step state
    pk.Nc = max(pk.Nc, I1.N);   pk.Nt = min(pk.Nt, I1.N);
    pk.M = max(pk.M, I1.M);
    pk.sig = max(pk.sig, abs(I1.N) ./ P.A3 + I1.M .* P.cI3);
    pk.Fn = max(pk.Fn, I1.Fn);  pk.Ft = max(pk.Ft, vecnorm(I1.Ft));
    pk.Fn_total = max(pk.Fn_total, sum(I1.Fn));
    pk.acc_g = max(pk.acc_g, norm(I1.acc_b) / 9.80665);
    pk.acc_axial_g = max(pk.acc_axial_g, I1.acc_b(1) / 9.80665);
    pk.tilt_deg = max(pk.tilt_deg, I1.tilt);
    a_lp = a_lp + alpha * (I1.acc_b - a_lp);
    fn_lp = fn_lp + alpha * (I1.Fn - fn_lp);
    pk.acc_lp_g = max(pk.acc_lp_g, norm(a_lp) / 9.80665);
    pk.Fn_lp = max(pk.Fn_lp, fn_lp);  pk.Fn_total_lp = max(pk.Fn_total_lp, sum(fn_lp));
    if opts.Log && mod(n, opts.LogEvery) == 0
        iLog = iLog + 1;
        H.t(iLog) = t;  H.z(iLog) = x(3);  H.vz(iLog) = x(6);  H.tilt(iLog) = I1.tilt;
        H.acc_g(iLog) = norm(I1.acc_b) / 9.80665;  H.acc_lp_g(iLog) = norm(a_lp) / 9.80665;
        H.N(:, :, iLog) = I1.N;  H.Fn(:, iLog) = I1.Fn;  H.Ft(:, iLog) = vecnorm(I1.Ft);
        H.s(:, :, iLog) = Z.s;  H.sp(:, iLog) = Z.sp;  H.dp(:, iLog) = Z.dp;
    end
    if n == nStep, break; end
    k2 = rates(x + dt/2 * k1, Z, t + dt/2, P);
    k3 = rates(x + dt/2 * k2, Z, t + dt/2, P);
    k4 = rates(x + dt * k3, Z, t + dt, P);
    x = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    x(7:10) = x(7:10) / norm(x(7:10));
    t = t + dt;
    [Z, dE] = plastic(x, Z, P);
    pk.E_core = pk.E_core + dE(1);  pk.E_pad = pk.E_pad + dE(2);
    if x(3) < -2 || pk.tilt_deg > 80, break; end   % fell through / over
    if mod(n, 500) == 0                           % settled: slow for SettleTime
        if t > 0.2 && norm(x(4:6)) < 0.02 && norm(x(11:13)) < 0.02
            settled = settled + 500 * dt;
            if settled >= opts.SettleTime, break; end
        else
            settled = 0;
        end
    end
end

pk.stroke_core = max(Z.s, [], 'all');          % main-strut core crush (m)
pk.stroke_pad = max(Z.sp);
pk.dent = max(Z.dp);
pk.buckling_ratio = pk.Nc ./ P.Pcr;            % peak compression / Euler load
pk.tipped = pk.tilt_deg > G.tip_angle_deg;
pk.E_in = E0;
R.peak = pk;
R.final.t = t;  R.final.x = x;  R.final.Z = Z;  R.final.E = energy(x, P);
if opts.Log
    keep = 1:iLog;
    for fn = string(fieldnames(H)).'
        a = H.(fn);
        if ndims(a) == 3, H.(fn) = a(:, :, keep); else, H.(fn) = a(:, keep); end
    end
    R.hist = H;
end
R.case = C;  R.surface = S.name;  R.core = G.core;  R.joints = G.joints;
end

% =====================================================================
function P = setup(G, S, C)
cg = [C.cg_x; 0; 0];
P.node_b = G.node - cg;
P.pad_b = G.pad - cg;
P.off_b = G.pad - G.node;                     % pad bottom relative to the foot node
P.up_b = G.upper - cg;                        % 3x3x4
P.L0 = G.L0;  P.k = G.k;  P.c = G.c;
P.A3 = repmat(G.A.', 1, 4);  P.cI3 = repmat((G.c_out ./ G.I).', 1, 4);
P.Pcr = G.Pcr;
P.Kb = G.Kb;  P.fixed = strcmp(G.joints, 'fixed');  P.Mmap = G.Mmap;
P.mf = G.m_foot;  P.c_foot = G.c_foot;
P.mb = C.mass - 4 * G.m_foot;
P.I = C.I;  P.Iinv = inv(C.I);
P.g = [0; 0; -9.80665];
P.T0 = C.thrust;  P.tau = C.thrust_tau;
P.prof = [];  P.t0 = 0;
if isfield(C, 'thrust_profile') && ~isempty(C.thrust_profile)
    P.prof = C.thrust_profile;  P.t0 = C.thrust_t0;
end
core = G.core;
P.coreMain = strcmp(core.placement, 'main');
P.corePad = strcmp(core.placement, 'pad');
P.cr = core;
ser = @(a, b) 1 ./ (1 ./ a + 1 ./ b);
if P.coreMain && isfield(core, 'k_elastic')
    P.k(1, :) = ser(P.k(1, :), core.k_elastic);       % core spring in series with the main strut
end
if P.corePad && isfield(core, 'k_elastic')
    S.k = ser(S.k, core.k_elastic);                   % core spring in series with the ground
    S.kbb = ser(S.kbb, core.k_elastic);
    S.Ku = ser(S.Ku, core.k_elastic);
end
P.S = S;
P.cg_n = 2 * S.zeta * sqrt(max(S.Ku, 1) * G.m_foot);   % ground damping
P.ct = 2 * 0.3 * sqrt(S.kt * G.m_foot);                  % friction spring damping
end

function F = coreForce(s, cr)
F = cr.plateau + (cr.peak - cr.plateau) .* max(0, 1 - s / cr.trigger) ...
    + cr.k_dens .* max(0, s - cr.s_max);
end

function F = soilForce(x, dp, S)
% elastic-plastic normal force at penetration x (>= 0), dent dp
if strcmp(S.type, 'soil')
    F = min(S.Ku .* (x - dp), S.kbb .* x);
else
    F = S.k .* x;
end
F = max(F, 0);
end

function [xd, I] = rates(x, Z, t, P)
p = x(1:3);  v = x(4:6);  q = x(7:10);  w = x(11:13);
f = reshape(x(14:25), 3, 4);  fv = reshape(x(26:37), 3, 4);
Rb = q2dcm(q);

rn = Rb.' * (f - p);                          % foot nodes relative to CG, body
vn = Rb.' * (fv - v) - cross(repmat(w, 1, 4), rn);   % relative velocity, body

% ---- struts ----
N = zeros(3, 4);  Fleg = zeros(3, 4);
for i = 1:3
    d = rn - squeeze(P.up_b(:, i, :));
    L = vecnorm(d);  u = d ./ L;
    del = P.L0(i, :) - L;                     % shortening
    ddot = -sum(u .* vn, 1);
    Ni = P.k(i, :) .* (del - Z.s(i, :)) + P.c(i, :) .* ddot;
    if i == 1 && P.coreMain
        Ni = min(Ni, coreForce(Z.s(1, :), P.cr));
    end
    N(i, :) = Ni;
    Fleg = Fleg + Ni .* u;                    % on the foot, body axes
end
Fleg = Fleg - P.c_foot * vn;
M = zeros(3, 4);
if P.fixed
    D = rn - P.node_b;
    for k = 1:4
        Fleg(:, k) = Fleg(:, k) - P.Kb(:, :, k) * D(:, k);
        for i = 1:3
            M(i, k) = max(norm(P.Mmap(:, :, 1, i, k) * D(:, k)), norm(P.Mmap(:, :, 2, i, k) * D(:, k)));
        end
    end
end

% ---- ground ----
S = P.S;
cp = f + Rb * P.off_b;                        % contact points
pen = -cp(3, :);
pr = -fv(3, :);                               % penetration rate
xs = pen - Z.sp;                              % ground-side compression
Fn = soilForce(xs, Z.dp, S);
on = xs > 0 & Fn > 0;
Fn = Fn + on .* P.cg_n .* pr;
Fn = max(Fn, 0) .* on;
if P.corePad
    Fn = min(Fn, coreForce(Z.sp, P.cr));
end
e = cp(1:2, :) - Z.a;
Ft = -S.kt * e - P.ct * fv(1:2, :);
mag = vecnorm(Ft);  lim = S.mu * Fn;
sc = min(1, lim ./ max(mag, eps));
Ft = Ft .* sc .* on;
Fg = [Ft; Fn];

% ---- accelerations ----
FlegI = Rb * Fleg;
fa = (FlegI + Fg) / P.mf + P.g;
if isempty(P.prof)
    T = P.T0 * exp(-t / P.tau);
else
    tc = t + P.t0;                            % time since the cutoff command
    if tc <= P.prof(1, 1), f = P.prof(2, 1);
    elseif tc >= P.prof(1, end), f = P.prof(2, end);
    else, f = interp1(P.prof(1, :), P.prof(2, :), tc);
    end
    T = P.T0 * f;
end
Fb = -sum(FlegI, 2) + T * Rb(:, 1);
va = Fb / P.mb + P.g;
Mb = sum(cross(rn, -Fleg), 2);
wd = P.Iinv * (Mb - cross(w, P.I * w));
qd = 0.5 * [-q(2:4).' * w; q(1) * w + cross(q(2:4), w)];
xd = [v; va; qd; wd; fv(:); fa(:)];

if nargout > 1
    I.N = N;  I.M = M;  I.Fn = Fn;  I.Ft = Ft;
    I.acc_b = Rb.' * (Fb / P.mb);             % specific force, body (what an IMU reads)
    I.tilt = acosd(min(1, Rb(3, 1)));
end
end

function [Z, dE] = plastic(x, Z, P)
p = x(1:3);  q = x(7:10);
f = reshape(x(14:25), 3, 4);
Rb = q2dcm(q);
rn = Rb.' * (f - p);
dE = [0 0];
if P.coreMain
    d = rn - squeeze(P.up_b(:, 1, :));
    del = P.L0(1, :) - vecnorm(d);
    s0 = Z.s(1, :);
    Fc = coreForce(s0, P.cr);
    k = P.k(1, :);
    yield = k .* (del - s0) > Fc;
    sn = del - Fc ./ k;
    dens = sn > P.cr.s_max;
    sn(dens) = (k(dens) .* del(dens) - P.cr.plateau + P.cr.k_dens * P.cr.s_max) ./ (k(dens) + P.cr.k_dens);
    sn = max(s0, sn);
    Z.s(1, yield) = sn(yield);
    dE(1) = sum(Fc .* (Z.s(1, :) - s0));
end
S = P.S;
cp = f + Rb * P.off_b;
pen = -cp(3, :);
if P.corePad
    sp0 = Z.sp;
    Fc = coreForce(sp0, P.cr);
    xs = pen - sp0;
    yield = soilForce(xs, Z.dp, S) > Fc;
    if strcmp(S.type, 'soil')
        xneed = max(Z.dp + Fc / S.Ku, Fc / S.kbb);
    else
        xneed = Fc / S.k;
    end
    Z.sp(yield) = max(sp0(yield), pen(yield) - xneed(yield));
    dE(2) = sum(Fc .* (Z.sp - sp0));
end
if strcmp(S.type, 'soil')
    xs = pen - Z.sp;
    over = S.Ku * (xs - Z.dp) > S.kbb * xs & xs > 0;
    Z.dp(over) = xs(over) - S.kbb * xs(over) / S.Ku;
end
% friction anchors
xs = pen - Z.sp;
Fn = soilForce(xs, Z.dp, S) .* (xs > 0);
c2 = cp(1:2, :);
e = c2 - Z.a;
Ftr = -S.kt * e;
mag = vecnorm(Ftr);  lim = S.mu * Fn;
slip = mag > lim;
Fs = Ftr ./ max(mag, eps) .* lim;
Z.a(:, slip) = c2(:, slip) + Fs(:, slip) / S.kt;
end

function E = energy(x, P)
v = x(4:6);  w = x(11:13);  fv = reshape(x(26:37), 3, 4);  f = reshape(x(14:25), 3, 4);
E = 0.5 * P.mb * (v.' * v) + 0.5 * w.' * P.I * w + 0.5 * P.mf * sum(fv.^2, 'all') ...
    - P.mb * P.g(3) * x(3) - P.mf * P.g(3) * sum(f(3, :));
end

function R = q2dcm(q)
a = q(1); b = q(2); c = q(3); d = q(4);
R = [a^2+b^2-c^2-d^2, 2*(b*c-a*d),     2*(b*d+a*c);
     2*(b*c+a*d),     a^2-b^2+c^2-d^2, 2*(c*d-a*b);
     2*(b*d-a*c),     2*(c*d+a*b),     a^2-b^2-c^2+d^2];
end

function q = dcm2q(R)
tr = trace(R);
if tr > 0
    s = 2 * sqrt(tr + 1);
    q = [s/4; (R(3,2)-R(2,3))/s; (R(1,3)-R(3,1))/s; (R(2,1)-R(1,2))/s];
else
    [~, i] = max(diag(R));
    switch i
        case 1
            s = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
            q = [(R(3,2)-R(2,3))/s; s/4; (R(1,2)+R(2,1))/s; (R(1,3)+R(3,1))/s];
        case 2
            s = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
            q = [(R(1,3)-R(3,1))/s; (R(1,2)+R(2,1))/s; s/4; (R(2,3)+R(3,2))/s];
        otherwise
            s = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
            q = [(R(2,1)-R(1,2))/s; (R(1,3)+R(3,1))/s; (R(2,3)+R(3,2))/s; s/4];
    end
end
q = q / norm(q);
end

function R = axang(a, th)
K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
end

function PLAN = gfold_onboard(r_ned, v_ned, m_now, IN, opts)
% GFOLD_ONBOARD  Fuel-optimal powered descent, solved onboard from the
%                current vehicle state.
%
%   PLAN = gfold_onboard(r_ned, v_ned, m_now, IN)
%   PLAN = gfold_onboard(r_ned, v_ned, m_now, IN, opts)
%
% Solves the descent leg (apex -> IN.mission.landing_alt, at rest) as a
% second-order cone program, using the same formulation as
% sims/GFOLD/gfold_descent_v3.m:
%
%   minimise  sum(Tm)*dt                      (total impulse ~ fuel)
%   s.t.      trapezoidal translational dynamics, mass frozen at m_now
%             ||T(:,k)|| <= Tm(k)              (SOC thrust magnitude)
%             Tmin <= Tm(k) <= Tmax            (throttle bounds)
%             T_up(k) >= Tm(k)*cos(theta_max)  (pointing / gimbal limit)
%             |vz(k)| <= vz_max                (descent rate limit)
%             alt(k) >= 0                      (no ground penetration)
%             sum(Tm)*dt <= m_prop*Isp*g0      (propellant budget)
%
% Every limit is read from IN (mission_inputs) -- nothing is hardcoded.
%
% INPUTS
%   r_ned  3x1  current position in NED  [m]   (z is DOWN, so alt = -z)
%   v_ned  3x1  current velocity in NED  [m/s]
%   m_now  1x1  current vehicle mass     [kg]
%   IN          struct from mission_inputs()
%   opts        optional struct:
%                 .N        number of intervals            (default 60)
%                 .tf_list  candidate flight times [s]     (default sweep;
%                           the first FEASIBLE entry in list order wins)
%                 .Isp      specific impulse [s]           (default 225)
%                 .verbose  print solver progress          (default false)
%
% OUTPUT  PLAN struct:
%   .ok      true if a feasible trajectory was found
%   .status  solver status string
%   .t       1xM relative time grid, t(1) = 0            [s]
%   .r_ned   3xM reference position, NED                 [m]
%   .v_ned   3xM reference velocity, NED                 [m/s]
%   .T_ned   3xM reference thrust vector, NED            [N]
%   .m       1xM reference mass                          [kg]
%   .tf      flight time of the returned plan            [s]
%   .fuel    propellant consumed by the returned plan    [kg]
%
% The caller converts PLAN into the LQR's xref / u_nom tables; see
% gfold_onboard_step.m.
%
% Solver: CVX when it is installed, otherwise coneprog (Optimization
% Toolbox R2020b+). Both paths solve the identical program.

if nargin < 5 || isempty(opts); opts = struct(); end

N       = getfielddef(opts, 'N', 60);
Isp     = getfielddef(opts, 'Isp', 225);
verbose = getfielddef(opts, 'verbose', false);

%% ---- limits, all from IN -------------------------------------------------
g0        = IN.const.g0;
T_nominal = IN.propulsion.nominal_thrust;
Tmin      = IN.propulsion.throttle_range(1) * T_nominal;
Tmax      = IN.propulsion.throttle_range(2) * T_nominal;
theta_max = IN.tvc.max_gimbal_angle;
vz_max    = IN.mission.max_descent_vel;
alt_f     = IN.mission.landing_alt;
m_dry     = IN.mass.dry_total;

m_prop = max(m_now - m_dry, 0);        % propellant still on board [kg]

%% ---- NED (z down) -> planner frame (z up) -------------------------------
r_ned = r_ned(:);  v_ned = v_ned(:);

r0 = [r_ned(1); r_ned(2); -r_ned(3)];
v0 = [v_ned(1); v_ned(2); -v_ned(3)];
rf = [0; 0; alt_f];
vf = [0; 0; 0];

g_vec = [0; 0; -g0];

%% ---- candidate flight times --------------------------------------------
% tf drives feasibility, and the binding limit is the Tmin THROTTLE FLOOR,
% not Tmax. The engine can never shut off, so the fastest the vehicle may
% descend is at a_down = g0 - Tmin/m. With 300 lbf nominal that floor is
% 533.8 N, which already exceeds vehicle weight below 54.4 kg -- at and
% below that mass no descent is feasible at any tf and this function
% correctly returns ok = false (FDIR F-23).
%
% So the ladder is scaled by the Tmin-limited fall time rather than the
% free-fall time, and is walked in ascending order. Fuel grows monotonically
% with tf inside the feasible band, so the first feasible tf is also the
% cheapest and the sweep stops there.
if isfield(opts, 'tf_list') && ~isempty(opts.tf_list)
    tf_list = opts.tf_list(:)';
else
    h0 = max(r0(3) - alt_f, 0);                      % height still to lose
    a_down = g0 - Tmin/m_now;                        % best downward accel
    if a_down <= 1e-3
        % Tmin alone holds or lifts the vehicle: it cannot descend.
        PLAN = emptyplan();
        PLAN.status = sprintf( ...
            'infeasible: Tmin=%.0f N >= weight=%.0f N at m=%.2f kg', ...
            Tmin, m_now*g0, m_now);
        return
    end
    t_ref   = sqrt(2*max(h0,1e-6)/a_down);
    tf_list = unique(round(t_ref * ...
        [0.9 1.0 1.1 1.25 1.4 1.6 1.8 2.1 2.5 3.0], 2));
    tf_list = tf_list(tf_list >= 0.5);
end

%% ---- solve ---------------------------------------------------------------
PLAN = emptyplan();

use_cvx = (exist('cvx_begin', 'file') == 2) || (exist('cvx_begin', 'file') == 6);

for tf = tf_list
    dt = tf / N;
    t  = linspace(0, tf, N+1);

    if use_cvx
        [r, v, T, Tm, status] = solve_cvx(N, dt, r0, v0, rf, vf, g_vec, ...
            m_now, Tmin, Tmax, theta_max, vz_max, m_prop, Isp, g0);
    else
        [r, v, T, Tm, status] = solve_coneprog(N, dt, r0, v0, rf, vf, g_vec, ...
            m_now, Tmin, Tmax, theta_max, vz_max, m_prop, Isp, g0);
    end

    if ~status.ok
        if verbose
            fprintf('[gfold_onboard] tf=%.2fs: %s\n', tf, status.msg);
        end
        continue
    end

    % mass reconstructed from the impulse profile (as in gfold_descent_v3)
    m_traj = zeros(1, N+1);
    m_traj(1) = m_now;
    for k = 1:N
        m_traj(k+1) = m_traj(k) - dt/(2*Isp*g0) * (Tm(k) + Tm(k+1));
    end

    PLAN.ok     = true;
    PLAN.status = status.msg;
    PLAN.t      = t;
    PLAN.r_ned  = [ r(1,:);  r(2,:); -r(3,:)];
    PLAN.v_ned  = [ v(1,:);  v(2,:); -v(3,:)];
    PLAN.T_ned  = [ T(1,:);  T(2,:); -T(3,:)];
    PLAN.m      = m_traj;
    PLAN.tf     = tf;
    PLAN.fuel   = m_now - m_traj(end);

    if verbose
        fprintf('[gfold_onboard] tf=%.2fs solved, fuel %.2f kg\n', tf, PLAN.fuel);
    end

    break   % ascending ladder -> first feasible tf is the cheapest
end

if verbose && ~PLAN.ok
    fprintf('[gfold_onboard] no feasible tf in [%.2f %.2f] s\n', ...
        min(tf_list), max(tf_list));
end

end

% ========================================================================
function [r, v, T, Tm, status] = solve_cvx(N, dt, r0, v0, rf, vf, g_vec, ...
    m0, Tmin, Tmax, theta_max, vz_max, m_prop, Isp, g0)

r = []; v = []; T = []; Tm = [];

cvx_begin quiet
    variable rr(3, N+1)
    variable vv(3, N+1)
    variable TT(3, N+1)
    variable TTm(1, N+1)

    minimize( sum(TTm) * dt )

    subject to
        rr(:,1) == r0;      vv(:,1) == v0;
        rr(:,N+1) == rf;    vv(:,N+1) == vf;

        for k = 1:N
            rr(:,k+1) == rr(:,k) + dt/2*(vv(:,k) + vv(:,k+1));
            vv(:,k+1) == vv(:,k) + dt/2*( ...
                (TT(:,k)/m0 + g_vec) + (TT(:,k+1)/m0 + g_vec));
        end

        for k = 1:N+1
            norm(TT(:,k), 2) <= TTm(k);
            TT(3,k) >= TTm(k) * cos(theta_max);
        end

        TTm >= Tmin;
        TTm <= Tmax;

        vv(3,:) >= -vz_max;
        vv(3,:) <=  vz_max;

        rr(3,:) >= 0;

        sum(TTm) * dt <= m_prop * Isp * g0;
cvx_end

status.ok  = strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved');
status.msg = cvx_status;

if status.ok
    r = rr; v = vv; T = TT; Tm = TTm;
end

end

% ========================================================================
function [r, v, T, Tm, status] = solve_coneprog(N, dt, r0, v0, rf, vf, g_vec, ...
    m0, Tmin, Tmax, theta_max, vz_max, m_prop, Isp, g0)
% Same program, assembled by hand for coneprog.
%
% Variable layout, 10 unknowns per node k = 1..Np:
%   x(off+1:off+3) = r(:,k)
%   x(off+4:off+6) = v(:,k)
%   x(off+7:off+9) = T(:,k)
%   x(off+10)      = Tm(k)          with off = (k-1)*10

r = []; v = []; T = []; Tm = [];

Np = N + 1;
n  = 10 * Np;

ir = @(k) (k-1)*10 + (1:3);
iv = @(k) (k-1)*10 + (4:6);
iT = @(k) (k-1)*10 + (7:9);
im = @(k) (k-1)*10 + 10;

%% objective: minimise dt * sum(Tm)
f = zeros(n, 1);
for k = 1:Np
    f(im(k)) = dt;
end

%% equality constraints
nEq = 12 + 6*N;
Aeq = spalloc(nEq, n, 12 + 30*N);
beq = zeros(nEq, 1);
row = 0;

% boundary conditions
Aeq(row+(1:3), ir(1))  = eye(3);  beq(row+(1:3)) = r0;  row = row + 3;
Aeq(row+(1:3), iv(1))  = eye(3);  beq(row+(1:3)) = v0;  row = row + 3;
Aeq(row+(1:3), ir(Np)) = eye(3);  beq(row+(1:3)) = rf;  row = row + 3;
Aeq(row+(1:3), iv(Np)) = eye(3);  beq(row+(1:3)) = vf;  row = row + 3;

for k = 1:N
    % r(k+1) - r(k) - dt/2*(v(k) + v(k+1)) = 0
    Aeq(row+(1:3), ir(k+1)) =  eye(3);
    Aeq(row+(1:3), ir(k))   = -eye(3);
    Aeq(row+(1:3), iv(k))   = -dt/2 * eye(3);
    Aeq(row+(1:3), iv(k+1)) = -dt/2 * eye(3);
    beq(row+(1:3))          = zeros(3,1);
    row = row + 3;

    % v(k+1) - v(k) - dt/(2*m0)*(T(k) + T(k+1)) = dt*g
    Aeq(row+(1:3), iv(k+1)) =  eye(3);
    Aeq(row+(1:3), iv(k))   = -eye(3);
    Aeq(row+(1:3), iT(k))   = -dt/(2*m0) * eye(3);
    Aeq(row+(1:3), iT(k+1)) = -dt/(2*m0) * eye(3);
    beq(row+(1:3))          = dt * g_vec;
    row = row + 3;
end

%% inequality constraints  (A*x <= b)
nIq = 4*Np + 1;
A   = spalloc(nIq, n, 4*Np + Np);
b   = zeros(nIq, 1);
row = 0;

ct = cos(theta_max);
for k = 1:Np
    jr = ir(k);  jv = iv(k);  jT = iT(k);  jm = im(k);

    % pointing: Tm*cos(theta_max) - T_up <= 0
    row = row + 1;
    A(row, jm)    =  ct;
    A(row, jT(3)) = -1;
    b(row)        =  0;

    % vz <= vz_max
    row = row + 1;
    A(row, jv(3)) =  1;
    b(row)        =  vz_max;

    % -vz <= vz_max
    row = row + 1;
    A(row, jv(3)) = -1;
    b(row)        =  vz_max;

    % no ground penetration: -alt <= 0
    row = row + 1;
    A(row, jr(3)) = -1;
    b(row)        =  0;
end

% propellant budget: dt*sum(Tm) <= m_prop*Isp*g0
row = row + 1;
for k = 1:Np
    A(row, im(k)) = dt;
end
b(row) = m_prop * Isp * g0;

%% bounds: Tmin <= Tm <= Tmax, everything else free
lb = -inf(n, 1);
ub =  inf(n, 1);
for k = 1:Np
    lb(im(k)) = Tmin;
    ub(im(k)) = Tmax;
end

%% second-order cones:  ||T(:,k)|| <= Tm(k)
% coneprog form: norm(Asc*x - bsc) <= d'*x - gamma
for k = Np:-1:1
    jT  = iT(k);
    Asc = sparse(1:3, jT, ones(1,3), 3, n);
    d   = sparse(im(k), 1, 1, n, 1);
    socs(k) = secondordercone(Asc, zeros(3,1), d, 0);
end

o = optimoptions('coneprog', 'Display', 'off');
[x, ~, exitflag] = coneprog(f, socs, A, b, Aeq, beq, lb, ub, o);

status.ok  = (exitflag == 1) && ~isempty(x);
status.msg = sprintf('coneprog exitflag %d', exitflag);

if ~status.ok
    return
end

r  = zeros(3, Np);
v  = zeros(3, Np);
T  = zeros(3, Np);
Tm = zeros(1, Np);
for k = 1:Np
    r(:,k)  = x(ir(k));
    v(:,k)  = x(iv(k));
    T(:,k)  = x(iT(k));
    Tm(k)   = x(im(k));
end

end

% ========================================================================
function PLAN = emptyplan()
PLAN = struct('ok', false, 'status', 'not solved', 't', [], ...
    'r_ned', [], 'v_ned', [], 'T_ned', [], 'm', [], 'tf', NaN, 'fuel', NaN);
end

function val = getfielddef(s, name, dflt)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = dflt;
end
end

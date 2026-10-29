function [best, results] = sweepBaffleDesign( ...
    prop_mass_max, prop_density, nu, tank_r, tank_h, opts)
% Sweep ring baffle count, width and thickness for one tank and return the
% best design meeting the minimum lateral and axial damping targets.
% opts.objective: 'count' (fewest baffles, then lightest; default) or
%                 'mass'  (lightest, then fewest baffles)
%
% Targets are enforced over fill fractions >= opts.min_fill (damping goes
% to the smooth-wall value as the tank empties for every design).
% Designs are screened on an n_mass grid; feasible ones are then
% re-checked in objective order on a finer n_verify grid, since axial
% damping dips sharply between baffles.
% If no design is feasible, the design with the smallest worst-case
% shortfall is returned with best.feasible = false.

if nargin < 6, opts = struct(); end
opts = setDefault(opts, 'zeta_lat_min', 0.04);
opts = setDefault(opts, 'zeta_ax_min',  0.30);
opts = setDefault(opts, 'min_fill',     0.10);
opts = setDefault(opts, 'Nb_vec',       1:40);
opts = setDefault(opts, 'w_frac_vec',   0.05:0.05:0.5);      % w / tank_r
opts = setDefault(opts, 't_vec',        [0.5 1 1.5 2 3]*1e-3); % m
opts = setDefault(opts, 'slosh_amp',    0.1 * tank_r);
opts = setDefault(opts, 'baffle_rho',   2700);                % Al6061
opts = setDefault(opts, 'n_mass',       300);
opts = setDefault(opts, 'n_verify',     3000);
opts = setDefault(opts, 'objective',    'count');

switch opts.objective
    case 'count', rank_keys = {'Nb', 'mass'};
    case 'mass',  rank_keys = {'mass', 'Nb'};
    otherwise, error('sweepBaffleDesign: objective must be ''count'' or ''mass''.');
end

mass_profile  = linspace(opts.min_fill * prop_mass_max, prop_mass_max, opts.n_mass);
verify_profile = linspace(opts.min_fill * prop_mass_max, prop_mass_max, opts.n_verify);

% thick/narrow candidates fall outside the TN D-3716 data; don't spam
warn_state = warning('off', 'baffleThicknessFactor:extrap');
restore_warn = onCleanup(@() warning(warn_state));

nN = numel(opts.Nb_vec); nW = numel(opts.w_frac_vec); nT = numel(opts.t_vec);
n  = nN * nW * nT;

Nb = zeros(n,1); w = zeros(n,1); t = zeros(n,1);
zeta_lat = zeros(n,1); zeta_ax = zeros(n,1); mass = zeros(n,1);

k = 0;
for iN = 1:nN
    for iW = 1:nW
        for iT = 1:nT
            k = k + 1;
            Nb(k) = opts.Nb_vec(iN);
            w(k)  = opts.w_frac_vec(iW) * tank_r;
            t(k)  = opts.t_vec(iT);

            [zeta_lat(k), zeta_ax(k)] = minDamping(mass_profile, prop_density, ...
                nu, tank_r, tank_h, Nb(k), w(k), t(k), opts.slosh_amp);
            mass(k) = opts.baffle_rho * Nb(k) * t(k) * pi * (tank_r^2 - (tank_r - w(k))^2);
        end
    end
end

shortfall = worstShortfall(zeta_lat, zeta_ax, opts);
feasible  = shortfall <= 0;
verified  = false(n,1);

% fine-grid check of screened designs, best first
[~, order] = sortrows(table(Nb, mass), rank_keys);
for k = order(feasible(order))'
    [zeta_lat(k), zeta_ax(k)] = minDamping(verify_profile, prop_density, ...
        nu, tank_r, tank_h, Nb(k), w(k), t(k), opts.slosh_amp);
    shortfall(k) = worstShortfall(zeta_lat(k), zeta_ax(k), opts);
    feasible(k)  = shortfall(k) <= 0;
    verified(k)  = true;
    if feasible(k)
        break
    end
end

results = table(Nb, w, t, mass, zeta_lat, zeta_ax, shortfall, feasible, verified);

if any(feasible & verified)
    cand = sortrows(results(feasible & verified, :), rank_keys);
else
    cand = sortrows(results, {'shortfall', 'mass'});
end

best = table2struct(cand(1, :));

end

function [lat_min, ax_min] = minDamping(mass_profile, prop_density, nu, ...
    tank_r, tank_h, Nb, w, t, slosh_amp)
lat_min = min(lateralDampingCalc(mass_profile, prop_density, nu, tank_r, ...
    w, Nb, tank_h, t, slosh_amp));
ax_min  = min(axialDampingCalc(mass_profile, prop_density, nu, tank_r, ...
    w, Nb, tank_h, t, slosh_amp));
end

function s = worstShortfall(zeta_lat, zeta_ax, opts)
% worst-case fractional shortfall against the two targets (<= 0 -> feasible)
s = max(1 - zeta_lat / opts.zeta_lat_min, 1 - zeta_ax / opts.zeta_ax_min);
end

function s = setDefault(s, name, val)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = val;
end
end

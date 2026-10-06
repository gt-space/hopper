
ox_density = 1141;   % kg/m^3
fu_density = 800;    % kg/m^3

fu_nu = 1.3e-6;  % m^2/s Fu
ox_nu = 0.17e-6; % m^2/s OX

tank_r = TANKS.singular.radius;
% as-built tank lengths (override sizer TANKS.singular.*.h)
ox_tank_h = 29.18 * 0.0254;   % LOx, m
fu_tank_h = 18.91 * 0.0254;   % RP-1, m

% 'fixed': closed-loop verified design below (default)
% 'sweep': analytic sizing to the damping targets (sweepBaffleDesign)
% baffle_mode / baffle_export may be preset by a caller (sweepBaffleSim)
if ~exist('baffle_mode', 'var'),   baffle_mode   = 'fixed'; end
if ~exist('baffle_export', 'var'), baffle_export = true;    end

baffle_opts.zeta_lat_min = 0.04;
baffle_opts.zeta_ax_min  = 0.15;
baffle_opts.min_fill     = 0.10;
baffle_opts.slosh_amp    = 0.1 * tank_r;   % wall slosh amplitude, m
baffle_opts.objective    = 'count';        % 'count' (fewest baffles) or 'mass'
baffle_opts.baffle_rho   = 2700;           % Al6061

% Closed-loop design (hopper_6dof_NED_v2, 2026-10-05): 0 flight
% divergences in 70 runs at 0.75x / 1x / 2x damping. Rings sit from the
% tank floor up to 'span' (max liquid level), not the full tank.
% Landing transient (last 1 m): >10 deg pitch in ~7% of runs vs ~1.5% for
% the analytic 33/21 design.
% Validated with the baffle mount's vertical members included; they added
% <= 0.0013 lateral damping (rodDampingCalc), well inside the 0.75x margin,
% so they are left out here.
ox_fixed = struct('Nb', 17, 'w', 0.5 * tank_r, 't', 0.5e-3, 'span', 14.59 * 0.0254);
fu_fixed = struct('Nb', 4,  'w', 0.5 * tank_r, 't', 0.5e-3, 'span', 17.48 * 0.0254);
in2m = 0.0254;

ox_mass_profile = linspace(0.1, IN.propulsion.oxidizer_mass, 200);
fu_mass_profile = linspace(0.1, IN.propulsion.fuel_mass, 200);

switch baffle_mode
    case 'sweep'
        [ox_baffle, ox_baffle_sweep] = sweepBaffleDesign(IN.propulsion.oxidizer_mass, ...
            ox_density, ox_nu, tank_r, ox_tank_h, baffle_opts);
        [fu_baffle, fu_baffle_sweep] = sweepBaffleDesign(IN.propulsion.fuel_mass, ...
            fu_density, fu_nu, tank_r, fu_tank_h, baffle_opts);
        ox_baffle.span = ox_tank_h;
        fu_baffle.span = fu_tank_h;
    case 'fixed'
        ox_baffle = ox_fixed; fu_baffle = fu_fixed;
        ox_baffle_sweep = []; fu_baffle_sweep = [];
    otherwise
        error('generateBaffleLUT: baffle_mode must be ''fixed'' or ''sweep''.');
end

ox_baffle_number    = ox_baffle.Nb;
ox_baffle_width     = ox_baffle.w;
ox_baffle_thickness = ox_baffle.t;
fu_baffle_number    = fu_baffle.Nb;
fu_baffle_width     = fu_baffle.w;
fu_baffle_thickness = fu_baffle.t;

[ox_lateral_damping_ratios] = lateralDampingCalc(ox_mass_profile, ox_density, ox_nu, tank_r, ox_baffle_width, ox_baffle_number, ox_baffle.span, ox_baffle_thickness, baffle_opts.slosh_amp);

[fu_lateral_damping_ratios] = lateralDampingCalc(fu_mass_profile, fu_density, fu_nu, tank_r, fu_baffle_width, fu_baffle_number, fu_baffle.span, fu_baffle_thickness, baffle_opts.slosh_amp);

[ox_axial_damping_ratios] = axialDampingCalc(ox_mass_profile, ox_density, ox_nu, tank_r, ox_baffle_width, ox_baffle_number, ox_baffle.span, ox_baffle_thickness, baffle_opts.slosh_amp);

[fu_axial_damping_ratios] = axialDampingCalc(fu_mass_profile, fu_density, fu_nu, tank_r, fu_baffle_width, fu_baffle_number, fu_baffle.span, fu_baffle_thickness, baffle_opts.slosh_amp);

tank_names = ["Ox", "Fuel"];
designs = {ox_baffle, fu_baffle};
lut_lat = {ox_lateral_damping_ratios, fu_lateral_damping_ratios};
lut_ax  = {ox_axial_damping_ratios, fu_axial_damping_ratios};
lut_m   = {ox_mass_profile, fu_mass_profile};
fprintf('\n=== Baffles (%s design; targets: lateral >= %.3f, axial >= %.3f for fill >= %.0f%%) ===\n', ...
    baffle_mode, baffle_opts.zeta_lat_min, baffle_opts.zeta_ax_min, 100 * baffle_opts.min_fill);
for i = 1:2
    d = designs{i};
    sel = lut_m{i} >= baffle_opts.min_fill * max(lut_m{i});
    d.mass = baffle_opts.baffle_rho * d.Nb * d.t * pi * (tank_r^2 - (tank_r - d.w)^2);
    fprintf('%s: Nb = %d, w = %.1f mm, t = %.1f mm, stack height = %.2f in, mass = %.3f kg, min lat = %.4f, min ax = %.4f\n', ...
        tank_names(i), d.Nb, 1e3 * d.w, 1e3 * d.t, d.span / in2m, d.mass, min(lut_lat{i}(sel)), min(lut_ax{i}(sel)));
end

figure
plot(ox_mass_profile, ox_lateral_damping_ratios)
yline(baffle_opts.zeta_lat_min, '--')
xline(baffle_opts.min_fill * IN.propulsion.oxidizer_mass, ':')
xlabel("Ox mass")
ylabel("Lateral Damping Ratio")
title("Lateral Damping Ratio over Ox mass, Number of Baffles: " + ox_baffle_number + " , Baffle width (m): " + ox_baffle_width + " , Baffle thickness (m): " + ox_baffle_thickness)

figure
plot(fu_mass_profile, fu_lateral_damping_ratios)
yline(baffle_opts.zeta_lat_min, '--')
xline(baffle_opts.min_fill * IN.propulsion.fuel_mass, ':')
xlabel("Fuel mass")
ylabel("Lateral Damping Ratio")
title("Lateral Damping Ratio over Fuel mass, Number of Baffles: " + fu_baffle_number + " , Baffle width(m): " + fu_baffle_width + " , Baffle thickness (m): " + fu_baffle_thickness)

figure
plot(ox_mass_profile, ox_axial_damping_ratios)
yline(baffle_opts.zeta_ax_min, '--')
xline(baffle_opts.min_fill * IN.propulsion.oxidizer_mass, ':')
xlabel("Ox mass")
ylabel("Axial Damping Ratio")
title("Axial Damping Ratio over Ox mass, Number of Baffles: " + ox_baffle_number + " , Baffle width (m): " + ox_baffle_width + " , Baffle thickness (m): " + ox_baffle_thickness)

figure
plot(fu_mass_profile, fu_axial_damping_ratios)
yline(baffle_opts.zeta_ax_min, '--')
xline(baffle_opts.min_fill * IN.propulsion.fuel_mass, ':')
xlabel("Fuel mass")
ylabel("Axial Damping Ratio")
title("Axial Damping Ratio over Fuel mass, Number of Baffles: " + fu_baffle_number + " , Baffle width(m): " + fu_baffle_width + " , Baffle thickness (m): " + fu_baffle_thickness)

% Export LUT for sim_setup (baffle_lut.csv) and full design record (.mat)
if baffle_export
    baffle_lut = table(ox_mass_profile(:), ox_lateral_damping_ratios(:), ox_axial_damping_ratios(:), ...
        fu_mass_profile(:), fu_lateral_damping_ratios(:), fu_axial_damping_ratios(:), ...
        'VariableNames', {'ox_mass_profile', 'ox_lateral_damping_ratios', 'ox_axial_damping_ratios', ...
                          'fu_mass_profile', 'fu_lateral_damping_ratios', 'fu_axial_damping_ratios'});
    writetable(baffle_lut, 'baffle_lut.csv');

    save('baffle_lut.mat', 'baffle_lut', 'ox_baffle', 'fu_baffle', 'baffle_opts', 'baffle_mode', ...
        'ox_baffle_sweep', 'fu_baffle_sweep', 'tank_r', 'ox_tank_h', 'fu_tank_h');

    fprintf('Baffle LUT written to baffle_lut.csv / baffle_lut.mat\n');
end

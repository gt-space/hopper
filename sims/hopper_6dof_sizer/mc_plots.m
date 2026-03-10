function mc_plots(results_file)

if nargin < 1; results_file = 'mc_results.mat'; end

load(results_file, 'results');

n = length(results);

% Filter to successful runs only
ok = arrayfun(@(r) r.status.success, results);
results = results(ok);
n_ok = length(results);

fprintf('Plotting %d / %d successful runs\n', n_ok, n);

% =========================================================================
% Interpolate all timeseries onto a common time grid for statistics
% =========================================================================
t_max  = min(arrayfun(@(r) r.timeseries.t(end), results));
t_grid = linspace(0, t_max, 500)';

alt_mat   = zeros(n_ok, 500);
mass_mat  = zeros(n_ok, 500);
ox_mat    = zeros(n_ok, 500);
fu_mat    = zeros(n_ok, 500);

for i = 1:n_ok
    ts = results(i).timeseries;
    alt_mat(i,:)  = interp1(ts.t, ts.altitude, t_grid, 'linear', 'extrap');
    mass_mat(i,:) = interp1(ts.t, ts.tot_mass, t_grid, 'linear', 'extrap');
    ox_mat(i,:)   = interp1(ts.t, ts.ox_mass,  t_grid, 'linear', 'extrap');
    fu_mat(i,:)   = interp1(ts.t, ts.fu_mass,  t_grid, 'linear', 'extrap');
end

% Stats
alt_mean  = mean(alt_mat,  1);
alt_std   = std(alt_mat,   1);
mass_mean = mean(mass_mat, 1);
mass_std  = std(mass_mat,  1);
ox_mean   = mean(ox_mat,   1);
fu_mean   = mean(fu_mat,   1);

% Nominal (scenario closest to mean of all inputs)
max_alts = arrayfun(@(r) r.flight.max_altitude, results);
[~, nom_idx] = min(abs(max_alts - median(max_alts)));

% Colors
c_spread  = [0.6 0.8 1.0];   % light blue for spread
c_sigma   = [0.2 0.5 0.9];   % blue for 1/2/3 sigma bands
c_mean    = [1.0 0.4 0.0];   % orange for mean
c_nominal = [0.1 0.8 0.1];   % green for nominal

% =========================================================================
% Figure 1 — Altitude vs Time
% =========================================================================
figure('Name','Altitude vs Time','Position',[100 100 800 500]);
hold on; grid on;

for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.altitude, 'Color', 'b', 'LineWidth', 0.5);
end

% 1-sigma band
fill([t_grid; flipud(t_grid)], ...
     [alt_mean - alt_std, fliplr(alt_mean + alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% 2-sigma band
fill([t_grid; flipud(t_grid)], ...
     [alt_mean - 2*alt_std, fliplr(alt_mean + 2*alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.08, 'EdgeColor', 'none');

% 3-sigma band
fill([t_grid; flipud(t_grid)], ...
     [alt_mean - 3*alt_std, fliplr(alt_mean + 3*alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.05, 'EdgeColor', 'none');

% Nominal and mean
ts_nom = results(nom_idx).timeseries;
plot(ts_nom.t, ts_nom.altitude, 'Color', c_nominal, 'LineWidth', 2.0);
plot(t_grid, alt_mean, 'Color', c_mean, 'LineWidth', 2.5);

% Target altitude line
yline(51, '--k', '51m Target', 'LabelVerticalAlignment', 'bottom');

xlabel('Time (s)'); ylabel('Altitude (m)');
title(sprintf('Altitude vs Time — %d MC Runs', n_ok));
legend({'MC Runs','1\sigma','2\sigma','3\sigma','Nominal','Mean'}, ...
    'Location','northwest');

% =========================================================================
% Figure 2 — 3D Trajectory
% =========================================================================
figure('Name','3D Trajectory','Position',[150 150 800 600]);
hold on; grid on;

for i = 1:n_ok
    ts = results(i).timeseries;
    plot3(ts.y, ts.z, ts.x, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end

ts_nom = results(nom_idx).timeseries;
plot3(ts_nom.y, ts_nom.z, ts_nom.x, 'Color', c_nominal, 'LineWidth', 2.0);

xlabel('Y (m)'); ylabel('Z (m)'); zlabel('Altitude X (m)');
title(sprintf('3D Trajectory — %d MC Runs', n_ok));
legend({'MC Runs','Nominal'}, 'Location','best');
view(45, 25);

% =========================================================================
% Figure 3 — Propellant Mass vs Time
% =========================================================================
figure('Name','Propellant Mass vs Time','Position',[200 200 800 500]);

subplot(3,1,1); hold on; grid on;
for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.tot_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end
fill([t_grid; flipud(t_grid)], ...
     [mass_mean - mass_std, fliplr(mass_mean + mass_std)]', ...
     c_sigma, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_grid, mass_mean, 'Color', c_mean, 'LineWidth', 2.5);
plot(results(nom_idx).timeseries.t, results(nom_idx).timeseries.tot_mass, ...
    'Color', c_nominal, 'LineWidth', 2.0);
ylabel('Total Prop (kg)'); title('Propellant Mass vs Time');
legend({'MC','1\sigma','Mean','Nominal'}, 'Location','northeast');

subplot(3,1,2); hold on; grid on;
for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.ox_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end
plot(t_grid, ox_mean, 'Color', c_mean, 'LineWidth', 2.5);
plot(results(nom_idx).timeseries.t, results(nom_idx).timeseries.ox_mass, ...
    'Color', c_nominal, 'LineWidth', 2.0);
ylabel('Ox Mass (kg)');

subplot(3,1,3); hold on; grid on;
for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.fu_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end
plot(t_grid, fu_mean, 'Color', c_mean, 'LineWidth', 2.5);
plot(results(nom_idx).timeseries.t, results(nom_idx).timeseries.fu_mass, ...
    'Color', c_nominal, 'LineWidth', 2.0);
ylabel('Fuel Mass (kg)'); xlabel('Time (s)');

% =========================================================================
% Figure 4 — Statistics Summary
% =========================================================================
max_alts   = arrayfun(@(r) r.flight.max_altitude,    results);
land_vels  = arrayfun(@(r) r.flight.landing_vel,     results);
dry_masses = arrayfun(@(r) r.vehicle.dry_mass,        results);
twrs       = arrayfun(@(r) r.vehicle.twr_initial,     results);

figure('Name','MC Statistics','Position',[250 250 900 600]);

subplot(2,2,1);
histogram(max_alts, 20, 'FaceColor', c_sigma);
xline(mean(max_alts),   'Color', c_mean,    'LineWidth', 2, 'Label', 'Mean');
xline(median(max_alts), 'Color', c_nominal, 'LineWidth', 2, 'Label', 'Median');
xline(51, '--k', 'LineWidth', 1.5, 'Label', 'Target');
xlabel('Max Altitude (m)'); ylabel('Count');
title(sprintf('Max Altitude  \\mu=%.1fm  \\sigma=%.1fm', mean(max_alts), std(max_alts)));

subplot(2,2,2);
histogram(land_vels, 20, 'FaceColor', c_sigma);
xline(mean(land_vels),   'Color', c_mean,    'LineWidth', 2, 'Label', 'Mean');
xline(median(land_vels), 'Color', c_nominal, 'LineWidth', 2, 'Label', 'Median');
xline(0.5, '--k', 'LineWidth', 1.5, 'Label', 'Limit');
xlabel('Landing Velocity (m/s)'); ylabel('Count');
title(sprintf('Landing Velocity  \\mu=%.2f  \\sigma=%.2f', mean(land_vels), std(land_vels)));

subplot(2,2,3);
histogram(dry_masses, 20, 'FaceColor', c_sigma);
xline(mean(dry_masses),   'Color', c_mean,    'LineWidth', 2, 'Label', 'Mean');
xline(median(dry_masses), 'Color', c_nominal, 'LineWidth', 2, 'Label', 'Median');
xlabel('Dry Mass (kg)'); ylabel('Count');
title(sprintf('Dry Mass  \\mu=%.1fkg  \\sigma=%.1fkg', mean(dry_masses), std(dry_masses)));

subplot(2,2,4);
histogram(twrs, 20, 'FaceColor', c_sigma);
xline(mean(twrs),   'Color', c_mean,    'LineWidth', 2, 'Label', 'Mean');
xline(median(twrs), 'Color', c_nominal, 'LineWidth', 2, 'Label', 'Median');
xline(1.0, '--k', 'LineWidth', 1.5, 'Label', 'Min TWR');
xlabel('Initial TWR'); ylabel('Count');
title(sprintf('Initial TWR  \\mu=%.2f  \\sigma=%.2f', mean(twrs), std(twrs)));

% =========================================================================
% Print 3-sigma summary to console
% =========================================================================
fprintf('\n========== 3-SIGMA SUMMARY (%d runs) ==========\n', n_ok);
fprintf('%-25s %8s %8s %8s %8s\n', 'Metric', 'Mean', 'Std', '-3sig', '+3sig');
fprintf('%s\n', repmat('-', 1, 61));

metrics = {'Max Altitude (m)',   max_alts;  ...
           'Landing Vel (m/s)',  land_vels; ...
           'Dry Mass (kg)',       dry_masses; ...
           'Initial TWR',         twrs};

for k = 1:size(metrics,1)
    name = metrics{k,1};
    vals = metrics{k,2};
    mu   = mean(vals);
    sig  = std(vals);
    fprintf('%-25s %8.3f %8.3f %8.3f %8.3f\n', name, mu, sig, mu-3*sig, mu+3*sig);
end
fprintf('%s\n', repmat('=', 1, 61));

end
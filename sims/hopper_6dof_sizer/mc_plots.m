function mc_plots(results_file)

if nargin < 1; results_file = 'mc_results.mat'; end

load(results_file, 'results');
load(results_file, 'results', 'nominal');

% Replace nom_idx references with nominal.timeseries directly
% e.g. instead of:
%   nominal.timeseries.altitude
% use:
%   nominal.timeseries.altitude
n = length(results);

% Filter to successful runs only
ok = arrayfun(@(r) r.status.success, results);
results = results(ok);
n_ok = length(results);

fprintf('Plotting %d / %d successful runs\n', n_ok, n);

c_mc  = [0.30 0.60 1.00];
c_nom = [0.10 0.80 0.20];
c_mean= [1.00 0.20 0.20];

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
%max_alts = arrayfun(@(r) r.flight.max_altitude, results);
%[~, nom_idx] = min(abs(max_alts - median(max_alts)));

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

    if i == 1
        h_mc = plot(ts.t, ts.altitude, 'Color', 'b', 'LineWidth', 0.5);
    else
        plot(ts.t, ts.altitude, 'Color', 'b', 'LineWidth', 0.5);
    end
end

% 1-sigma band
h_sigma1 = fill([t_grid; flipud(t_grid)], ...
     [alt_mean - alt_std, fliplr(alt_mean + alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% 2-sigma band
h_sigma2 = fill([t_grid; flipud(t_grid)], ...
     [alt_mean - 2*alt_std, fliplr(alt_mean + 2*alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.08, 'EdgeColor', 'none');

% 3-sigma band
h_sigma3 = fill([t_grid; flipud(t_grid)], ...
     [alt_mean - 3*alt_std, fliplr(alt_mean + 3*alt_std)]', ...
     c_sigma, 'FaceAlpha', 0.05, 'EdgeColor', 'none');

% Nominal and mean
ts_nom = nominal.timeseries;
h_nom = plot(ts_nom.t, ts_nom.altitude, 'Color', c_nominal, 'LineWidth', 2.0);

h_mean = plot(t_grid, alt_mean, 'Color', c_mean, 'LineWidth', 2.5);

% Target altitude
yline(51, '--k', '51m Target', 'LabelVerticalAlignment', 'bottom');

xlabel('Time (s)');
ylabel('Altitude (m)');
title(sprintf('Altitude vs Time — %d MC Runs', n_ok));

legend([h_mc h_sigma1 h_sigma2 h_sigma3 h_nom h_mean], ...
       {'MC Runs','1\sigma','2\sigma','3\sigma','Nominal','Mean'}, ...
       'Location','northwest');
% =========================================================================
% Figure 2 — 3D Trajectory
% =========================================================================

figure('Name','3D Trajectory','Position',[150 150 800 600]);
hold on; grid on;

for i = 1:n_ok
    ts = results(i).timeseries;

    if i == 1
        h_mc = plot3(ts.x, ts.y, -ts.z, 'Color', 'b', 'LineWidth', 0.5);
    else
        plot3(ts.x, ts.y, -ts.z, 'Color', 'b', 'LineWidth', 0.5);
    end
end

ts_nom = nominal.timeseries;
h_nom = plot3(ts_nom.x, ts_nom.y, -ts_nom.z, 'Color', c_nominal, 'LineWidth', 2.0);

xlabel('Latitude (X)');
ylabel('Longitude (Y)');
zlabel('Altitude (Z)');
title(sprintf('3D Trajectory — %d MC Runs', n_ok));

legend([h_mc h_nom], {'MC Runs','Nominal'}, 'Location','best');

view(45, 25);
% =========================================================================
% Figure 3 — Propellant Mass vs Time
% =========================================================================
figure('Name','Propellant Mass vs Time','Position',[200 200 800 500]);

subplot(3,1,1); hold on; grid on;

for i = 1:n_ok
    ts = results(i).timeseries;

    if i == 1
        h_mc = plot(ts.t, ts.tot_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
    else
        plot(ts.t, ts.tot_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
    end
end

h_sigma = fill([t_grid; flipud(t_grid)], ...
     [mass_mean - mass_std, fliplr(mass_mean + mass_std)]', ...
     c_sigma, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

h_mean = plot(t_grid, mass_mean, 'Color', c_mean, 'LineWidth', 2.5);

h_nom = plot(nominal.timeseries.t, ...
             nominal.timeseries.tot_mass, ...
             'Color', c_nominal, 'LineWidth', 2.0);

ylabel('Total Prop (kg)');
title('Propellant Mass vs Time');

legend([h_mc h_sigma h_mean h_nom], ...
       {'MC Runs','1\sigma','Mean','Nominal'}, ...
       'Location','northeast');


subplot(3,1,2); hold on; grid on;
for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.ox_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end
plot(t_grid, ox_mean, 'Color', c_mean, 'LineWidth', 2.5);
plot(nominal.timeseries.t, nominal.timeseries.ox_mass, ...
    'Color', c_nominal, 'LineWidth', 2.0);
ylabel('Ox Mass (kg)');


subplot(3,1,3); hold on; grid on;
for i = 1:n_ok
    ts = results(i).timeseries;
    plot(ts.t, ts.fu_mass, 'Color', [c_spread 0.3], 'LineWidth', 0.5);
end
plot(t_grid, fu_mean, 'Color', c_mean, 'LineWidth', 2.5);
plot(nominal.timeseries.t, nominal.timeseries.fu_mass, ...
    'Color', c_nominal, 'LineWidth', 2.0);
ylabel('Fuel Mass (kg)');
xlabel('Time (s)');

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


% ============================
% Figure 5 - Thrust plot
% ============================

figure('Name','Thrust vs Time','Position',[250 50 900 450]);
hold on; grid on;
 
thrust_mat = zeros(n_ok, 500);
for i = 1:n_ok
    ts = results(i).timeseries;
    h_mc = plot(ts.t, ts.thrust, 'Color', [c_mc 0.4], 'LineWidth', 0.8);
    thrust_mat(i,:) = interp1(ts.t, ts.thrust, t_grid, 'linear', 'extrap');
end
h_mean = plot(t_grid, mean(thrust_mat,1), 'Color', c_mean, 'LineWidth', 3.0);
h_nom = plot(nominal.timeseries.t, nominal.timeseries.thrust, ...
    'Color', c_nom, 'LineWidth', 3.0);
 
xlabel('Time (s)','FontSize',12); 
ylabel('Thrust (N)','FontSize',12);
title(sprintf('Thrust vs Time — %d MC Runs', n_ok),'FontSize',14);
%legend({'MC Runs','Mean','Nominal'},'Location','best','FontSize',10);

legend([h_mc h_mean h_nom], {'MC Runs','Mean','Nominal'}, ...
       'Location','best','FontSize',10);


% ==================================
% Figure 6 - Velocity plot vs. Time
% ==================================
figure('Name','Velocity vs Time','Position',[300 100 900 600]);

vel_fields = {'vx','vy','vz_raw'};
vel_labels = {'Vx — North (m/s)','Vy — East (m/s)','Vz — Down (m/s)'};

for s = 1:3
    subplot(3,1,s); hold on; grid on;

    vel_mat = zeros(n_ok, 500);

    for i = 1:n_ok
        ts = results(i).timeseries;

        if i == 1
            h_mc = plot(ts.t, ts.(vel_fields{s}), 'Color', [c_mc 0.4], 'LineWidth', 0.8);
        else
            plot(ts.t, ts.(vel_fields{s}), 'Color', [c_mc 0.4], 'LineWidth', 0.8);
        end

        vel_mat(i,:) = interp1(ts.t, ts.(vel_fields{s}), t_grid, 'linear', 'extrap');
    end

    h_mean = plot(t_grid, mean(vel_mat,1), 'Color', c_mean, 'LineWidth', 3.0);

    h_nom = plot(nominal.timeseries.t, ...
                 nominal.timeseries.(vel_fields{s}), ...
                 'Color', c_nom, 'LineWidth', 3.0);

    ylabel(vel_labels{s},'FontSize',10);

    if s == 1
        legend([h_mc h_mean h_nom], {'MC Runs','Mean','Nominal'}, ...
               'Location','best','FontSize',9);
        title('Velocity vs Time','FontSize',13);
    end

    if s == 3
        xlabel('Time (s)','FontSize',12);
    end
end
% ==================================
% Figure 7 - Attitude plot vs. Time
% ==================================

figure('Name','Attitude vs Time','Position',[350 150 900 600]);

att_fields = {'roll','pitch','yaw'};
att_labels = {'Roll (deg)','Pitch (deg)','Yaw (deg)'};

for s = 1:3
    subplot(3,1,s); hold on; grid on;

    att_mat = zeros(n_ok, 500);

    for i = 1:n_ok
        ts = results(i).timeseries;
        data_deg = rad2deg(ts.(att_fields{s}));

        if i == 1
            h_mc = plot(ts.t, data_deg, 'Color', [c_mc 0.4], 'LineWidth', 0.8);
        else
            plot(ts.t, data_deg, 'Color', [c_mc 0.4], 'LineWidth', 0.8);
        end

        att_mat(i,:) = interp1(ts.t, data_deg, t_grid, 'linear', 'extrap');
    end

    h_mean = plot(t_grid, mean(att_mat,1), 'Color', c_mean, 'LineWidth', 3.0);

    h_nom = plot(nominal.timeseries.t, ...
        rad2deg(nominal.timeseries.(att_fields{s})), ...
        'Color', c_nom, 'LineWidth', 3.0);

    ylabel(att_labels{s},'FontSize',10);

    if s == 1
        legend([h_mc h_mean h_nom], {'MC Runs','Mean','Nominal'}, ...
               'Location','best','FontSize',9);
        title('Attitude vs Time','FontSize',13);
    end

    if s == 3
        xlabel('Time (s)','FontSize',12);
    end
end

% ================================
% Figure 8 - 2d trajectory with pad
% =================================
pad_radius = 0.5; % 0.5m diameter
land_x = arrayfun(@(r) r.timeseries.x(end), results);  % North
land_y = arrayfun(@(r) r.timeseries.y(end), results);  % East
on_pad = sqrt(land_x.^2 + land_y.^2) <= pad_radius;
figure('Name','Landing Scatter','Position',[400 100 650 650]);
hold on; grid on; axis equal;
% Pad circle
theta = linspace(0, 2*pi, 300);
fill(pad_radius*cos(theta), pad_radius*sin(theta), ...
    [0.85 0.85 0.85], 'EdgeColor','k', 'LineWidth', 2.0, 'FaceAlpha', 0.4);
% Landing points — color by on/off pad
h_off = scatter(land_x(~on_pad), land_y(~on_pad), 60, c_mean, 'filled', ...
    'MarkerEdgeColor','k');
h_on  = scatter(land_x(on_pad),  land_y(on_pad),  60, c_nom,  'filled', ...
    'MarkerEdgeColor','k');
% CEP circle (50th percentile radial error)
radial = sqrt(land_x.^2 + land_y.^2);
cep    = median(radial);
%plot(cep*cos(theta), cep*sin(theta), '--', 'Color', [0.3 0.6 1.0], 'LineWidth', 1.5);
% Mean landing point
plot(mean(land_x), mean(land_y), '+k', 'MarkerSize', 14, 'LineWidth', 2.5);
xlabel('North (m)','FontSize',12);
ylabel('East (m)','FontSize',12);
title(sprintf('Landing Scatter — %d MC Runs', n_ok),'FontSize',14);
legend([h_on h_off], ...
    {sprintf('On pad (%d)', sum(on_pad)), sprintf('Off pad (%d)', sum(~on_pad))}, ...
    'Location','best','FontSize',10);
% Pad and CEP annotations
text(0, pad_radius*1.1, sprintf('Pad \\phi=%.1fm', pad_radius*2), ...
    'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
%text(cep*0.7, cep*0.7, sprintf('CEP=%.2fm', cep), ...
    %'FontSize',9,'Color',[0.3 0.6 1.0]);
fprintf('\nLanding accuracy: %d / %d on pad (%.1f%%)\n', ...
    sum(on_pad), n_ok, 100*sum(on_pad)/n_ok);
%fprintf('CEP (median radial error): %.3f m\n', cep);
fprintf('Max radial error:          %.3f m\n', max(radial));

% new plots for OX MASS AND WIND

% =========================================================================
% Figure 9 — Oxidizer Mass vs Time
% =========================================================================
figure('Name','Oxidizer Mass vs Time','Position',[450 100 900 450]);
hold on; grid on;

ox_mat2 = zeros(n_ok, 500);
for i = 1:n_ok
    ts = results(i).timeseries;
    h_mc = plot(ts.t, ts.ox_mass, 'Color', [c_mc 0.4], 'LineWidth', 0.8);
    ox_mat2(i,:) = interp1(ts.t, ts.ox_mass, t_grid, 'linear', 'extrap');
end
h_mean = plot(t_grid, mean(ox_mat2,1), 'Color', c_mean, 'LineWidth', 3.0);
h_nom  = plot(nominal.timeseries.t, nominal.timeseries.ox_mass, ...
    'Color', c_nom, 'LineWidth', 3.0);

xlabel('Time (s)','FontSize',12); ylabel('Ox Mass (kg)','FontSize',12);
title(sprintf('Oxidizer Mass vs Time — %d MC Runs', n_ok),'FontSize',14);
legend([h_mc h_mean h_nom], {'MC Runs','Mean','Nominal'},'Location','best','FontSize',10);

% =========================================================================
% Figure 10 — Wind per MC Run (constant values as horizontal lines)
% =========================================================================
figure('Name','Wind MC Values','Position',[500 150 900 500]);

subplot(2,1,1); hold on; grid on;
for i = 1:n_ok
    uw = results(i).scenario.uwind;
    h_mc = plot([results(i).timeseries.t(1) results(i).timeseries.t(end)], ...
        [uw uw], 'Color', [c_mc 0.4], 'LineWidth', 1.5);
end
h_nom = plot([nominal.timeseries.t(1) nominal.timeseries.t(end)], ...
    [nominal.scenario.uwind nominal.scenario.uwind], 'Color', c_nom, 'LineWidth', 3.0);
ylabel('uwind (m/s)','FontSize',10);
title('Wind Values per MC Run','FontSize',13);
legend([h_mc h_nom], {'MC Runs','Nominal'},'Location','best','FontSize',9);

subplot(2,1,2); hold on; grid on;
for i = 1:n_ok
    vw = results(i).scenario.vwind;
    h_mc = plot([results(i).timeseries.t(1) results(i).timeseries.t(end)], ...
        [vw vw], 'Color', [c_mc 0.4], 'LineWidth', 1.5);
end
h_nom = plot([nominal.timeseries.t(1) nominal.timeseries.t(end)], ...
    [nominal.scenario.vwind nominal.scenario.vwind], 'Color', c_nom, 'LineWidth', 3.0);
ylabel('vwind (m/s)','FontSize',10);
xlabel('Time (s)','FontSize',12);
legend([h_mc h_nom], {'MC Runs','Nominal'},'Location','best','FontSize',9);

% new ox mass plot :)
% =========================================================================
% Figure 11 — Ox Mass Histogram (sampled input values)
% =========================================================================
figure('Name','Ox Mass Distribution','Position',[550 200 700 450]);
hold on; grid on;

ox_sampled = arrayfun(@(r) r.scenario.ox_mass, results);

histogram(ox_sampled, 15, 'FaceColor', c_mc, 'EdgeColor', 'w');
%xline(mean(ox_sampled),   'r-', 'LineWidth', 2.5, 'Label', 'Mean');
%xline(median(ox_sampled), 'g-', 'LineWidth', 2.5, 'Label', 'Median');
xline(nominal.scenario.ox_mass, 'k--', 'LineWidth', 2.0, 'Label', 'Nominal');

xlabel('Oxidizer Mass (kg)', 'FontSize', 12);
ylabel('Count', 'FontSize', 12);
title(sprintf('Ox Mass Distribution — %d MC Runs   \\mu=%.2f   \\sigma=%.2f', ...
    n_ok, mean(ox_sampled), std(ox_sampled)), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);


end
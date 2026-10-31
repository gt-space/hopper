% compare_sims  Compare reduced-physics configurations of hopper_6dof with the full one.
%
% Builds its own inputs (buildSimInput), so nothing needs to be in the workspace.
% Runs four configurations of the harness (they replace the old _simplified and
% _minimal model copies):
%   1) original          every axis at its default
%   2) no slosh          SLOSH=off
%   3) no slosh, calm    SLOSH=off, wind speeds 0
%   4) minimal           SLOSH=off, MASS_PROPS=constant (mass/inertia/CG frozen), calm
% then prints timings and the max difference for each logged signal, and plots overlays.
%
% READING THE TIMINGS: the configurations do not all run for the same amount of
% simulated time (the sim stops on ground contact), so total wall time is not a fair
% comparison. Use the 'us/sample' column, which is time per simulated step.
% There is also roughly 20% machine noise run to run, so ignore small differences.
%
% READING THE DIFFERENCES: roll/yaw differences near 180 or 360 are angle wrapping.
% The minimal configuration holds mass at its wet value, so it needs more thrust as
% the flight goes on and will drift from the original by design.

signals    = {'x', 'y', 'z', 'u', 'velocity', 'roll', 'pitch', 'yaw', 'thrust', 'cg', 'ox_mass', 'fuel_mass'};

nRuns = 3;   % runs per configuration; run 1 is warm-up (compile/cache) and is excluded from the average

W = hopperWorkspace();
Wcalm = W; Wcalm.IN.wind.uwind = 0; Wcalm.IN.wind.vwind = 0;
% name, variant choices, workspace
cfg = { ...
    'original',         variantConfig(),                                         W; ...
    'no slosh',         variantConfig('SLOSH', 'off'),                           W; ...
    'no slosh, calm',   variantConfig('SLOSH', 'off'),                           Wcalm; ...
    'minimal',          variantConfig('SLOSH', 'off', 'MASS_PROPS', 'constant'), Wcalm};

outs = cell(1, size(cfg, 1));
Ts   = cell(1, size(cfg, 1));
for k = 1:size(cfg, 1)
    in = buildSimInput('hopper_6dof', [], cfg{k,3}, cfg{k,2});
    fprintf('Running %s (%d runs)...\n', cfg{k,1}, nRuns);
    [outs{k}, Ts{k}] = timed_runs(in, nRuns);
end

report_perf(cfg, Ts, outs);

for k = 2:size(cfg, 1)
    report(sprintf('%s vs original', cfg{k,1}), outs{1}, outs{k}, signals);
end

styles = {'k-', 'b--', 'r-', 'g-.'};
for s = signals
    present = find(cellfun(@(o) has_sig(o, s{1}), outs));
    if numel(present) < 2, continue; end
    figure('Name', s{1}); hold on; grid on;
    for k = present
        plot_sig(outs{k}, s{1}, styles{min(k, numel(styles))}, cfg{k,1});
    end
    title(s{1}, 'Interpreter', 'none'); xlabel('time (s)'); legend('Location', 'best');
end

% ===================== helpers =====================
function [out, T] = timed_runs(in, nRuns)
% Run a SimulationInput nRuns times, recording wall-clock time, Simulink execution time and
% MATLAB memory growth. Returns the last run's output and a struct of timings.
    T.wall = zeros(1, nRuns); T.exec = nan(1, nRuns); T.mem = nan(1, nRuns);
    for r = 1:nRuns
        m0 = mem_used();
        tic;
        out = sim(in);
        T.wall(r) = toc;
        T.mem(r) = mem_used() - m0;                                  % MB, crude
        try
            T.exec(r) = out.SimulationMetadata.TimingInfo.ExecutionElapsedWallTime;
        catch
        end
    end
end

function mb = mem_used()
    try
        u = memory; mb = u.MemUsedMATLAB / 2^20;
    catch
        mb = nan;                                                    % memory() is Windows-only
    end
end

function report_perf(cfg, Ts, outs)
    fprintf('\n=== Performance (run 1 = warm-up, excluded from the average) ===\n');
    fprintf('%-18s %9s %9s %9s %9s %9s %10s %8s\n', ...
        'config', 'first(s)', 'avg(s)', 'exec avg', 'samples', 'sim end', 'us/sample', 'MB grown');
    base = nan;
    for k = 1:size(cfg, 1)
        T = Ts{k}; n = numel(T.wall);
        w = T.wall(min(2, n):end);
        e = T.exec(min(2, n):end);
        ns = numel(outs{k}.tout);
        us = 1e6 * mean(e, 'omitnan') / ns;                          % time per simulated step
        fprintf('%-18s %9.2f %9.2f %9.2f %9d %9.2f %10.1f %8.0f\n', cfg{k,1}, T.wall(1), mean(w), ...
            mean(e, 'omitnan'), ns, outs{k}.tout(end), us, mean(T.mem(min(2, n):end), 'omitnan'));
        if k == 1
            base = us;
        else
            fprintf('%-18s -> %.0f%% of the original per simulated step\n', '', 100 * us / base);
        end
    end
end

function tf = has_sig(out, s)
    tf = ~isempty(out.find(s));
end

function [t, d] = get_sig(out, s)
    v = out.get(s);
    if isa(v, 'timeseries'), t = v.Time; d = v.Data;
    elseif isstruct(v),      t = v.time; d = v.signals.values;
    else, error('Unsupported logged type for %s', s); end
    d = reshape(d, numel(t), []);
end

function plot_sig(out, s, style, name)
    [t, d] = get_sig(out, s);
    plot(t, d(:,1), style, 'DisplayName', [name ' (col 1)']);
end

function report(label, a, b, signals)
    fprintf('\n--- %s ---\n', label);
    for k = 1:numel(signals)
        s = signals{k};
        if ~has_sig(a, s) || ~has_sig(b, s)
            fprintf('%-10s (not logged in both)\n', s);
            continue;
        end
        [ta, da] = get_sig(a, s); [tb, db] = get_sig(b, s);
        if numel(ta) < 2 || numel(tb) < 2
            fprintf('%-10s (constant, too few samples to compare)\n', s);
            continue;
        end
        tt = max(min(ta), min(tb)) : median(diff(ta)) : min(max(ta), max(tb));
        ia = interp1(ta, da, tt(:), 'linear');
        ib = interp1(tb, db, tt(:), 'linear');
        fprintf('%-10s max|diff| = %-12.4g  (signal peak = %.4g, compared over %.2f s)\n', ...
            s, max(abs(ia(:) - ib(:))), max(abs(ia(:))), tt(end) - tt(1));
    end
end

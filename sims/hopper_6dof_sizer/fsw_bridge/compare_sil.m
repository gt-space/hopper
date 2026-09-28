function compare_sil(varargin)
%COMPARE_SIL  Compare Rust SIL runs against the Simulink closed-loop run.
%   compare_sil('build/sil_luna.csv', 'build/sil_simsign.csv', ...) compares
%   each Rust run (hopper_sil output, paths relative to fsw_bridge) with the
%   gust-free Simulink reference that verify_env_c saved. Run verify_env_c
%   first so the reference exists, and build hopper_sil against the same
%   gust-free code (HOPPER_ENV_CODE = build/verify_nogust/hopper_env_ert_rtw).

here = fileparts(mfilename('fullpath'));
ref  = load(fullfile(here, 'build', 'verify_nogust', 'verify_env_c.mat'), 't', 'U', 'Xs', 'zSL');
dt   = ref.t(2) - ref.t(1);
fprintf('Simulink reference: %d steps, ends t = %.3f s, max altitude %.3f m\n', ...
    numel(ref.t), ref.t(end), -min(ref.zSL));

groups = {'position (m)', 1:3; 'velocity (m/s)', 4:6; 'body rate (rad/s)', 7:9; 'quaternion', 10:13};
for f = 1:numel(varargin)
    R = readmatrix(fullfile(here, varargin{f}), 'NumHeaderLines', 1);
    k = round(R(:, 1) / dt) + 1;
    keep = k <= numel(ref.t);
    R = R(keep, :); k = k(keep);
    X = R(:, 2:14); U = R(:, 15:18); z = R(:, 20);

    fprintf('\n== %s: ends t = %.3f s, max altitude %.3f m ==\n', varargin{f}, R(end, 1), -min(z));
    fprintf('%-18s %12s %12s\n', 'signal', 'max|diff|', 'at t (s)');
    for g = 1:size(groups, 1)
        [m, j] = max(max(abs(X(:, groups{g, 2}) - ref.Xs(k, groups{g, 2})), [], 2));
        fprintf('%-18s %12.4g %12.3f\n', groups{g, 1}, m, R(j, 1));
    end
    unames = {'thrust cmd (N)', 'TVC pitch cmd', 'TVC yaw cmd', 'RCS cmd'};
    for i = 1:4
        [m, j] = max(abs(U(:, i) - ref.U(k, i)));
        fprintf('%-18s %12.4g %12.3f\n', unames{i}, m, R(j, 1));
    end
    D = max(abs(X - ref.Xs(k, :)), [], 2);
    fprintf('first time state difference exceeds: ');
    for tol = [1e-12 1e-9 1e-6 1e-3]
        j = find(D > tol, 1);
        if isempty(j), fprintf('%g never; ', tol); else, fprintf('%g @ %.3f s; ', tol, R(j, 1)); end
    end
    fprintf('\n');
end
end

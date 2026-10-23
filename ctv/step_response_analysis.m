% Run sim_setup_ctv

% Time Ascent --> Early landing
t_main = [7, 12, 17, 25];

% Close to and within unstable region
t_land = [26.0, 26.3, 26.45, 26.65, 26.75, 27.0];

pairs = [6 1; 4 2; 8 2; 5 3; 9 3; 7 4];
pair_names = {'Vz \leftarrow Thrust', 'Vx \leftarrow \delta_p', 'Q \leftarrow \delta_p', ...
    'Vy \leftarrow \delta_y', 'R \leftarrow \delta_y', 'P \leftarrow RCS'};

plotStepSet(t_main, pairs, pair_names, tgrid, A, B, K1, 'Main Mission', 25);
plotStepSet(t_land, pairs, pair_names, tgrid, A, B, K1, 'Landing', 4);

function plotStepSet(t_points, pairs, pair_names, tgrid, A, B, K1, phase_label, t_end)
for p = 1:size(pairs,1)
    out_idx = pairs(p,1);
    in_idx  = pairs(p,2);
    figure;
    theme(gcf, "light");
    hold on;
    for k = 1:numel(t_points)
        A_k = lookupMat(tgrid, A, t_points(k));
        B_k = lookupMat(tgrid, B, t_points(k));
        K1_k = lookupMat(tgrid, K1, t_points(k));
        A_cl_k = A_k - B_k*K1_k; 
        sys_k = ss(A_cl_k, B_k, eye(13), zeros(13,4));
        step(sys_k(out_idx, in_idx), t_end); 
    end
    hold off;
    legend(arrayfun(@(t) sprintf('%.2fs', t), t_points, 'UniformOutput', false));
    title([phase_label ': ' pair_names{p}]);
end
end

function [T, nZero] = modeTable(A_cl)
    e = eig(A_cl);
    e = e(imag(e) >= 0);                   % one of each conjugate pair

    tol = 1e-9;
    sigma = real(e);  sigma(abs(sigma) < tol) = 0;
    omega = imag(e);  omega(abs(omega) < tol) = 0;

    zeroMode = (sigma == 0) & (omega == 0);  % numerical-noise / marginal modes
    nZero = sum(zeroMode);
    sigma = sigma(~zeroMode);
    omega = omega(~zeroMode);

    wn     = sqrt(sigma.^2 + omega.^2);
    zeta   = -sigma ./ wn;
    period = 2*pi ./ omega;
    t_half = log(2) ./ abs(sigma);

    eigenvalues = strings(numel(sigma), 1);
    for i = 1:numel(sigma)
        if omega(i) == 0
            eigenvalues(i) = sprintf('%.4g', sigma(i));
        else
            eigenvalues(i) = sprintf('%.4g ± j%.4g', sigma(i), omega(i));
        end
    end

    T = table(eigenvalues, wn, zeta, omega, period, t_half, ...
        'VariableNames', {'eigenvalues','omega_n','zeta','omega_d','T_period','t_half'});
end

t_tables = [12, 25, 27];

for k = 1:numel(t_tables)
    [~, idx] = min(abs(tgrid - t_tables(k)));
    A_cl_k = A(:,:,idx) - B(:,:,idx)*K1(:,:,idx);

    [T, nZero] = modeTable(A_cl_k);
    fprintf('\nt = %.2f s (snapshot at %.3f s): %d mode(s) at zero\n', ...
        t_tables(k), tgrid(idx), nZero);
    disp(T)

    exportModeTable(T, sprintf('mode_table_%gs.csv', t_tables(k)));
end

function exportModeTable(T, filename)
n = height(T);

header = ["Mode Name", "Eigenvalue(s) σ ± jω", "Nat. freq. ωn (rad/s)", ...
    "Damping ratio ζ", "Damped freq. ωd (rad/s)", "Time period T (s)", ...
    "Time-to-half / double (s)"];

fid = fopen(filename, 'w', 'n', 'UTF-8');
fwrite(fid, [239 187 191], 'uint8');
fprintf(fid, '%s\n', strjoin('"' + header + '"', ','));

for i = 1:n
    row = ["", T.eigenvalues(i), fmtNum(T.omega_n(i)), fmtNum(T.zeta(i)), ...
        fmtNum(T.omega_d(i)), fmtNum(T.T_period(i)), fmtNum(T.t_half(i))];
    fprintf(fid, '%s\n', strjoin('"' + row + '"', ','));
end
fclose(fid);
end

function s = fmtNum(x)
if isinf(x) || isnan(x)
    s = "n/a";
else
    s = string(sprintf('%.4g', x));
end
end
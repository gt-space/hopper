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
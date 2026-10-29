function T = compareSloshML(slosh_ml_dir, design)
% Compare the hopper_6dof_NED_v2 slosh model against SLOSH-ML
% (https://github.com/lgst-lab/slosh_ml, v1.0.0) for the as-built LOx and
% RP-1 tanks.
%
%   slosh_ml_dir - path to a clone of lgst-lab/slosh_ml
%   design       - optional struct with fields ox, fu (Nb, w, t, span);
%                  defaults to the generateBaffleLUT fixed design
%
% Lateral mode: frequency, slosh mass and pendulum length from the model's
% Slosh/MATLAB Function equations vs SLOSH-ML's Lomen/Dodge solver.
% Damping: lateralDampingCalc vs the SLOSH-ML correlations (Stephens,
% Mikishev smooth wall; Bauer, Coogan ring baffles at eta = 0.1 R).

here = fileparts(mfilename('fullpath'));
addpath(fullfile(slosh_ml_dir, 'source', 'v1.0.0', 'functions'));
addpath(fileparts(fileparts(here)));   % hopper_6dof_sizer

g    = 9.81;
in2m = 0.0254;
R    = 4.25 * in2m;
if nargin < 2
    design.ox = struct('Nb', 17, 'w', 0.5 * R, 't', 0.5e-3, 'span', 14.59 * in2m);
    design.fu = struct('Nb', 4,  'w', 0.5 * R, 't', 0.5e-3, 'span', 17.48 * in2m);
end
tanks = [struct('name', "LOx",  'rho', 1141, 'nu', 0.17e-6, 'H', 29.18 * in2m, 'mmax', 16, 'd', design.ox), ...
         struct('name', "RP-1", 'rho', 800,  'nu', 1.3e-6,  'H', 18.91 * in2m, 'mmax', 13, 'd', design.fu)];

xi = 1.841;
rows = {};
for T_ = tanks
    d = T_.d;
    for frac = [0.05 0.10 0.25 0.50 0.75 1.00]
        m = frac * T_.mmax;
        h = m / T_.rho / (pi * R^2);
        s = sloshMLHeadless(R, T_.H, h, T_.rho, g);

        % hopper_6dof_NED_v2 Slosh/MATLAB Function
        w_model = sqrt(xi * g / R * tanh(xi * h / R));
        m_model = m * (2 * R * tanh(xi * h / R)) / (xi * (xi^2 - 1) * h);
        l_model = R / (xi * tanh(xi * h / R));

        z_hop   = lateralDampingCalc(m, T_.rho, T_.nu, R, d.w, d.Nb, d.span, d.t, 0.1 * R);
        z_steph = damping_ratio_cylinder("Cylindrical (Stephens)", h, R, g, T_.nu);
        z_mik   = damping_ratio_cylinder("Cylindrical (Mikishev, NASA SP-106)", h, R, g, T_.nu);
        zb  = linspace(0, d.span, d.Nb);
        sub = zb(zb < h);
        if isempty(sub)
            z_bauer = 0; z_coog = 0;
        else
            spacing = d.span / max(d.Nb - 1, 1);
            top_d   = h - max(sub);
            % SLOSH-ML returns zeta / (eta/R)^0.5
            z_bauer = sqrt(0.1) * damping_ratio_cylinder_baffles("Cylindrical w/ Baffles (Bauer)", ...
                h, R, g, T_.nu, d.w, 1.55, spacing, numel(sub), top_d);
            z_coog  = sqrt(0.1) * damping_ratio_cylinder_baffles("Cylindrical w/ Baffles (Coogan)", ...
                h, R, g, T_.nu, d.w, 1.55, spacing, numel(sub), top_d);
        end

        rows(end+1, :) = {T_.name, frac, h / R, s.omega1, w_model, s.m1, m_model, s.L1, l_model, ...
            s.H_hinge1, s.H_spring1, s.m0, z_hop, z_steph, z_mik, z_bauer, z_coog}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'tank', 'fill', 'h_R', 'w_SLOSHML', 'w_model', ...
    'm1_SLOSHML', 'm1_model', 'L1_SLOSHML', 'L1_model', 'Hhinge_SLOSHML', 'Hbob_SLOSHML', 'm0_SLOSHML', ...
    'zeta_model', 'zeta_Stephens', 'zeta_Mikishev', 'zeta_Bauer', 'zeta_Coogan'});
T.w_err_pct  = 100 * (T.w_model  - T.w_SLOSHML)  ./ T.w_SLOSHML;
T.m1_err_pct = 100 * (T.m1_model - T.m1_SLOSHML) ./ T.m1_SLOSHML;
T.L1_err_pct = 100 * (T.L1_model - T.L1_SLOSHML) ./ T.L1_SLOSHML;
end

function [xref, unom, gfold_flag] = gfold_replan(t, x, m, ...
                                        xref_stored, unom_stored, tgrid)
%#codegen
% GFOLD_REPLAN  Body of the "GFOLD Replan" MATLAB Function block in
%               ctv/hopper_6dof_NED_v3_ctv.slx.
%
% This file is the single source of truth for that block's script:
% ctv/gnc/add_onboard_gfold.m reads it with fileread() and installs it, so
% edit it here and re-run add_onboard_gfold.m rather than editing the block.
%
% The real work is done by gfold_onboard_step, which is declared extrinsic
% so it runs interpreted and may call CVX / coneprog and hold persistent
% state. Outputs are preassigned so Simulink can infer their size and type
% without being able to see inside the extrinsic call -- xref and unom are
% seeded from the stored tables, which fixes them at N x 13 and N x 4.
%
% INPUTS
%   t             sim time                                      [s]
%   x             full 13-state [x y z vx vy vz P Q R q0..q3], NED
%   m             current vehicle mass                          [kg]
%   xref_stored   precomputed reference table   (Target_Trajectory_full)
%   unom_stored   precomputed nominal control table  (unom)
%   tgrid         reference time grid                (tgrid)
%
% OUTPUTS
%   xref, unom    tables to hand to Gain_Selector, with the descent portion
%                 replaced by the onboard G-FOLD solution
%   gfold_flag    0 = nominal, 1/2 = FDIR F-23 (see gfold_onboard_step)

coder.extrinsic('gfold_onboard_step');

% seed outputs so their size/type is known at compile time
xref       = xref_stored;
unom       = unom_stored;
gfold_flag = 0;

[xref, unom, gfold_flag] = gfold_onboard_step(t, x, m, ...
                                xref_stored, unom_stored, tgrid);

end

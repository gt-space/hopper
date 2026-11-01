function K = fdir_catalog_cols()
% FDIR_CATALOG_COLS  Column indices of the numeric fault table.
%
%#codegen

K.N          = 1;
K.M          = 2;
K.CUSUM      = 3;
K.FILTER     = 4;    % 1 N_OF_M, 2 CONSECUTIVE, 3 CUSUM
K.CLEAR      = 5;    % 1 NEVER_IN_FLIGHT, 2 AUTO_CLEAR, 3 GROUND_ONLY
K.PHASE_MASK = 6;    % bit k-1 set: active in phase k (fdir_modes)
K.RESPONSE   = 7;    % 0..4 R0..R4, 5 NO_ARM, 6 R4_ABOVE_HCUT, 7 STUCK_THRUST
K.ESC_TO     = 8;    % -1 none, else 0..4
K.ESC_AFTER  = 9;    % s latched before escalating
K.INHIBITED  = 10;   % 1: not evaluated (e.g. needs NAVIGATION)
K.TARGET     = 11;   % fdir_targets index
K.NCOL       = 11;

K.F_N_OF_M = 1; K.F_CONSECUTIVE = 2; K.F_CUSUM = 3;
K.C_NEVER_IN_FLIGHT = 1; K.C_AUTO_CLEAR = 2; K.C_GROUND_ONLY = 3;
K.R_NO_ARM = 5; K.R_R4_ABOVE_HCUT = 6; K.R_STUCK_THRUST = 7;
end

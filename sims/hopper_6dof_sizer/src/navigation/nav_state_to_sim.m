function [x, dcm_be] = nav_state_to_sim(xn, cg, F)
%NAV_STATE_TO_SIM  The nav filter's 20-state as the controller's 13-state.
%   [x, dcm_be] = nav_state_to_sim(xn, cg, F)   inverse of nav_state_from_sim
%     xn   nav estimate (20): CG position and velocity, nav body rates,
%          nav quaternion, biases
%     cg   Mass_Properties CG [x; y; z], sim body (m)
%     F    nav_frames
%   x       [pos NED; vel NED; rates (sim body); quaternion (sim)] (13), in
%           the harness's convention (sim_pos_offset): what the controller
%           and FDIR read in place of truth
%   dcm_be  NED -> sim body from the estimate
%
%#codegen

xn = xn(:);
q = xn(10:13);
q = q / max(norm(q), eps);
dcm_be = F.R' * nav_dcm_from_quat(q);
x = zeros(13, 1);
x(1:3) = xn(1:3) - [0; 0; F.heightCOM] + sim_pos_offset(cg);
x(4:6) = xn(4:6);
x(7:9) = F.R' * xn(7:9);
x(10:13) = nav_quat_from_dcm(dcm_be);       % the sim's quat2dcm convention is the same A(q)
end

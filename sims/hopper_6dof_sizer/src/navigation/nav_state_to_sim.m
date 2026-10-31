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
% q_sim = q_nav * q_c (Hamilton), q_c the quaternion of the fixed frame
% change: the same attitude as dcm_be, but it keeps the filter's sign. A
% DCM round trip (w >= 0) flips the sign wherever w crosses zero (the CTV
% turning about its long axis), and nav_error and anything else comparing
% quaternions directly then see a full turn (the LQR scheduler aligns the
% sign itself)
qc = nav_quat_from_dcm(F.R');
x(10:13) = [q(1)*qc(1) - q(2:4)'*qc(2:4); q(1)*qc(2:4) + qc(1)*q(2:4) + cross(q(2:4), qc(2:4))];
end

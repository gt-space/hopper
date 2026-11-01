function xn = nav_state_from_sim(x, dcm_be, cg, F)
%NAV_STATE_FROM_SIM  Sim truth state as the nav filter's 20-state.
%   xn = nav_state_from_sim(x, dcm_be, cg, F)
%     x       the controller's 13-state [pos NED; vel NED; rates (sim
%             body); quaternion (sim)]: CG position shifted to the base
%             (sim_pos_offset) and CG velocity
%     dcm_be  NED -> sim body
%     cg      Mass_Properties CG [x; y; z], sim body (m)
%     F       nav_frames
%   xn      [pos NED (CG, from its pad rest point); vel NED (CG);
%           rates (nav body); quaternion (nav, scalar first);
%           accel bias; gyro bias; baro bias] (20), biases zero
%
%   For the filter's initial state and for checking the filter: a perfect
%   estimate is xn.
%
%#codegen

x = x(:);
xn = zeros(20, 1);
xn(1:3) = x(1:3) - sim_pos_offset(cg) + [0; 0; F.heightCOM];
xn(4:6) = x(4:6);
xn(7:9) = F.R * x(7:9);
xn(10:13) = nav_quat_from_dcm(F.R * dcm_be);   % NED -> nav body
end

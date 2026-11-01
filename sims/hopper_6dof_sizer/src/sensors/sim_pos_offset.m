function d = sim_pos_offset(cg)
%SIM_POS_OFFSET  What the harness adds to the CG position to make x(1:3).
%   d = sim_pos_offset(cg)   cg: Mass_Properties CG [x; y; z], sim body (m)
%
%   The controller's state x(1:3) is the 6DOF block's CG position (NED)
%   plus d: the harness's Sum10 / Sum12 / Sum13 blocks subtract the CG
%   offset as if the vehicle stood on the pad (body x up, z north, y east),
%   without rotating it by the attitude. x(4:6) is the CG velocity as is.
%   The nav adapter undoes exactly this, so a perfect estimate maps back
%   onto the sim's own state.
%
%#codegen

d = [-cg(3); -cg(2); cg(1)];
end

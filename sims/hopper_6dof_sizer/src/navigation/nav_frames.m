function F = nav_frames(dcm_pad, cg_pad)
%NAV_FRAMES  Fixed relation between the sim's frames and the nav filter's.
%   F = nav_frames(dcm_pad, cg_pad)
%     dcm_pad  sim DCM (NED to sim body) of the vehicle standing on the pad
%     cg_pad   CG height above the base on the pad (m, sim body x): cg_init
%
%   F.R          sim body -> nav body (3x3): v_nav = F.R * v_sim
%   F.heightCOM  nav position origin above the pad (m): the CG's resting
%                height, so the nav position is 0 on the pad and the ground
%                plane sits at nav z = +heightCOM (lidar_update_model)
%
%   Sim body: x up the thrust axis; positions in the vehicle config are
%   from the base reference point.
%   Nav body (subsystems/navigation): z down the long axis, origin at the
%   CG; position is the CG's, NED, from its resting point on the pad. The
%   nav attitude is the identity on the pad (their hover tests), so
%   F.R = dcm_pad': on the pad the nav body axes are N, E, D, and sim +x
%   (up) is nav -z.
%
%#codegen

F.R = dcm_pad';
F.heightCOM = cg_pad;
end

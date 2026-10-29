function P = navParams(W, S, dt)
%NAVPARAMS  Frames, sigmas and sensor geometry the nav adapter uses.
%   P = navParams(W, S, dt)   W from hopperWorkspace (cg_init), S the
%                             Sensors block's sim_sensor_params, dt the
%                             filter step
%   P.F (nav_frames), P.S, P.geom (lidar beams and mounts, GPS antennas:
%   sim body, from the base), P.dt
%
%   navWorkspace builds the filter's parameters from it, and FDIR's
%   estimator checks (fdir_nav_monitors) read the same, so both use the
%   same measurement geometry.

dcm_pad = angle2dcm(0, pi/2, 0);     % the 6DOF block's eul_0 = [0 pi/2 0]: body x up
P.F = nav_frames(dcm_pad, W.cg_init);
P.S = S;
[P.geom.lidar_dir, P.geom.lidar_mount] = lidar_geometry(S.lidar);
P.geom.gps_antenna = S.gps.antenna;
P.dt = dt;
end

%% sim_setup_ctv
% Populates the base workspace for hopper_6dof_NED_v2_ctv.slx.
%
% v2 IS THE ACTIVE CTV MODEL. It is the only ctv model with evidence of ever
% having been built (ctv/hopper_6dof_NED_v2_ctv.slxc plus
% ctv/slprj/.../hopper_6dof_NED_v2_ctv/), and its cached variable list was
% used to verify that this script produces everything it needs.
%
% v2 reads Target_Trajectory and its Gain_Selector does
% reshape(Kflat, 13, 4), i.e. the FULL 13-state controller with a 52-column
% K1flat. That is what LinerizationMaster_ctv builds, so that is what runs
% below.
%
% v3 (hopper_6dof_NED_v3_ctv.slx) wants the REDUCED 12-state controller
% instead: Target_Trajectory_full and a 48-column K1flat. LQTrack_Master_ctv
% builds that, and it is kept in the repo for whenever v3 is picked up, but
% it is NOT part of this script. v3 also still has an unresolved reference to
% IN.tanks.lateral_damping / IN.tanks.axial_damping, which v2 does not use.
%
% ONE-TIME MODEL CHANGES (run once each, they edit the .slx):
%   >> cd gnc
%   >> disable_roll_actuator          % required: see note below
%   >> add_onboard_gfold              % optional: onboard G-FOLD descent
%   >> cd ..
%
% disable_roll_actuator is REQUIRED. IN.rcs.max_thrust = 0 on its own makes
% the EDF Dynamics subsystem compute 0/0 = NaN, which propagates into the
% 6DOF moments. That script removes the subsystem structurally instead.

clearvars;

main_ctv

load_lookup_ctv

cg_moi_test_ctv

LinerizationMaster_ctv

generateBaffleLUT_ctv

%prop_system_testing;

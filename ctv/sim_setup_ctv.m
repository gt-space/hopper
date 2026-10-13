%% sim_setup_ctv
% Populates the base workspace for hopper_6dof_NED_v3_ctv.slx.
%
% v3 reads a Constant block whose value is Target_Trajectory_full and its
% Gain_Selector does reshape(Kflat, 12, 4)', i.e. it needs the REDUCED
% 12-state controller: a 48-column K1flat. LinerizationMaster_ctv.m builds
% the full 13-state controller (Target_Trajectory, 52-column K1flat), which
% is what hopper_6dof_NED_v2_ctv.slx needs instead. This script therefore
% runs LQTrack_Master_ctv for v3; run LinerizationMaster_ctv by hand if you
% want to drive the v2 model.
%
% KNOWN GAP (pre-existing, not fixed here): the v3 model also references
% IN.tanks.lateral_damping and IN.tanks.axial_damping in its Slosh
% subsystem, but ctv/inputs/mission_inputs.m leaves both commented out and
% generateBaffleLUT_ctv does not assign them either, so the model cannot
% resolve them. sims/hopper_6dof_sizer/inputs/mission_inputs.m has the same
% two lines uncommented and marked "MUST BE UNCOMMENTED". See the PR notes.

clearvars;

main_ctv

load_lookup_ctv

cg_moi_test_ctv

LQTrack_Master_ctv

generateBaffleLUT_ctv

%prop_system_testing;

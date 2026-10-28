function A = variantAxes()
%VARIANTAXES  The variant axes of models/hopper_6dof.slx and their choices.
%   A = variantAxes() returns a struct with one field per axis. Each field
%   lists the axis's choices; the first one is the default, which matches
%   hopper_6dof_NED_v2.
%
%   The model selects a choice with the variant control expression
%   VARIANT_<AXIS> == k, where k is the choice's position in this list, so
%   add new choices at the end (the build script and variantConfig both
%   read this list). See docs/VARIANTS.md.
%
%   Axis        Choices
%   VEHICLE     hopper        the hopper (inputs/mission_inputs.m)
%               ctv           controls test vehicle (inputs/ctv_inputs.m);
%                             a parameter set, not a model block. Defaults
%                             SLOSH to off and MASS_PROPS to constant.
%   SLOSH       full          slosh pendulums (subsystems/slosh/Slosh.slx)
%               off           no slosh forces, moments or CG shift
%   PROPULSION  lut           thrust command -> mdot lookup (Thrust_to_Mdot_LUT)
%               fluid_network #26's Simscape feed system and engine
%   GROUND      none          no ground reaction (v2)
%               strut         #27's strut spring/damper per leg
%   MASS_PROPS  varying       propellant drains: mass, CG, MOI and dI/dt change
%               constant      no propellant flow: mass, CG and MOI stay at
%                             their initial values, dI/dt = 0 (the engine
%                             gimbal still shifts the CG by < 1 mm)
%   CONTROLLER  lqr_sched     time-scheduled LQR along the reference (v2)
%               lqi           LQI-DEV's fixed-gain LQI hover regulator
%                             (experimental: toy plant, not tuned)
%               lqr_sched_legacy  the March lqr/LQ_6DOF_V3 scheduler, without
%                             the quaternion sign fix

A.VEHICLE    = ["hopper", "ctv"];
A.SLOSH      = ["full", "off"];
A.PROPULSION = ["lut", "fluid_network"];
A.GROUND     = ["none", "strut"];
A.MASS_PROPS = ["varying", "constant"];
A.CONTROLLER = ["lqr_sched", "lqi", "lqr_sched_legacy"];
end

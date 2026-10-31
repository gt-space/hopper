function parameterizeRollActuation()
%PARAMETERIZEROLLACTUATION  Roll disturbance and fan torque limit per vehicle.
%   Two hopper values reached every vehicle:
%   - models/hopper_6dof.slx's random roll disturbance (Random Gusts:
%     Bernoulli gate x Gaussian torque into the body roll moment) had its
%     variance typed in as 20 (N m)^2. It now reads
%     IN.vehicle.disturbance.roll_torque_var (vehicles/<vehicle>/..._config.m).
%   - subsystems/actuators/EDF_Dynamics.slx's only active limit, Saturation1
%     on the output torque, was +-IN.rcs.max_thrust: a fan thrust in N used
%     as a torque in N m (16 N m where 16 N on the 0.108 m arm gives
%     1.73 N m). It is now +-IN.rcs.max_thrust * IN.vehicle.rcs.arm. The
%     hopper's nominal flight peaks at 0.7 N m, so it is unchanged.
%   Rerunning is harmless.

m = 'hopper_6dof';
load_system(m);
rs = [m '/Random' newline 'Source'];
var = 'IN.vehicle.disturbance.roll_torque_var';
if strcmp(get_param(rs, 'VarVal'), var)
    fprintf('%s: roll disturbance already parameterized\n', m);
else
    assert(strcmp(get_param(rs, 'VarVal'), '20'), 'unexpected Random Source variance %s', get_param(rs, 'VarVal'));
    set_param(rs, 'VarVal', var);
    save_system(m);
    fprintf('%s: roll disturbance variance now %s\n', m, var);
end
close_system(m, 0);

f = 'EDF_Dynamics';
load_system(f);
sat = [f '/Saturation1'];
lim = 'IN.rcs.max_thrust * IN.vehicle.rcs.arm';
if strcmp(get_param(sat, 'UpperLimit'), lim)
    fprintf('%s: torque limit already max_thrust x arm\n', f);
else
    assert(strcmp(get_param(sat, 'UpperLimit'), 'IN.rcs.max_thrust'), 'unexpected limit %s', get_param(sat, 'UpperLimit'));
    set_param(sat, 'UpperLimit', lim, 'LowerLimit', ['-(' lim ')']);
    save_system(f);
    fprintf('%s: torque limit now +-%s\n', f, lim);
end
close_system(f, 0);
end

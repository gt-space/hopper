function G = scopeGroups()
%SCOPEGROUPS  The scope groups in the harness Logging block, and what each shows.
%   G = scopeGroups()   struct array: name, buses, displays
%
%   Each group is a Variant Subsystem Logging/Scopes_<name>, off unless
%   named in variantConfig('SCOPES', ...); a group that's off isn't
%   compiled, so it costs nothing. Turning one on opens its scope when the
%   simulation starts.
%
%   displays: one row per scope display, {title, signals}. A signal is
%   '<bus>.<element>' with an optional transform after a colon:
%     :deg      times 180/pi
%     :diag     a 3x3 matrix's diagonal
%     :i1-3     elements 1 to 3
%   or 'ref.<variable>', a workspace timeseries read with From Workspace
%   (references the setup already makes). Signals in one display are
%   overlaid.
%
%   To add or change a group: edit this list, then rerun
%   tools/refactor/buildLoggingScopes (it rebuilds every group).

G = struct('name', {}, 'buses', {}, 'displays', {});
G(end+1) = grp('trajectory', {
    'Position x, y (m)',                {'dyn_bus.x', 'dyn_bus.y'}
    'Altitude z vs reference (m)',      {'dyn_bus.z', 'ref.z_ref_ts'}
    'Vertical speed vs reference (m/s)', {'dyn_bus.w', 'ref.vz_ref_ts'}
    'Velocity NED (m/s)',               {'dyn_bus.vel_ned'}
    });
G(end+1) = grp('attitude', {
    'Roll, pitch, yaw (deg)',           {'dyn_bus.phi:deg', 'dyn_bus.theta:deg', 'dyn_bus.psi:deg'}
    'Body rates (deg/s)',               {'dyn_bus.w_body:deg'}
    'Quaternion',                       {'dyn_bus.quat'}
    });
G(end+1) = grp('propulsion', {
    'Thrust: actual, command, request, profile (N)', {'prop_bus.thrust', 'act_bus.thrust_cmd', 'act_bus.thrust_req', 'ref.T_ts'}
    'Propellant mass: ox, fuel (kg)',   {'prop_bus.ox_mass', 'prop_bus.fu_mass'}
    'Vehicle mass (kg)',                {'prop_bus.mass'}
    });
G(end+1) = grp('actuation', {
    'TVC pitch, yaw',                   {'act_bus.tvc_pitch', 'act_bus.tvc_yaw'}
    'Body force (N)',                   {'act_bus.force'}
    'Body moment (N m)',                {'act_bus.moment'}
    'Controller command u',             {'gnc_bus.u_cmd'}
    });
G(end+1) = grp('mass_props', {
    'CG (m)',                           {'prop_bus.cg'}
    'CG with slosh (m)',                {'prop_bus.cg_slosh'}
    'Engine CG (m)',                    {'prop_bus.engine_cg'}
    'MOI diagonal (kg m^2)',            {'prop_bus.moi:diag'}
    'dI/dt diagonal (kg m^2/s)',        {'prop_bus.moi_dot:diag'}
    });
G(end+1) = grp('slosh', {
    'Ox slug position (m)',             {'dyn_bus.ox_slosh_pos'}
    'Fuel slug position (m)',           {'dyn_bus.fu_slosh_pos'}
    'Slug masses: ox, fuel (kg)',       {'dyn_bus.ox_slosh_mass', 'dyn_bus.fu_slosh_mass'}
    'Specific force, body (m/s^2)',     {'dyn_bus.f_body'}
    'CG acceleration (m/s^2)',          {'dyn_bus.cg_accel'}
    });
G(end+1) = grp('environment', {
    'Wind force (N)',                   {'env_bus.wind_force'}
    'Wind moment (N m)',                {'env_bus.wind_moment'}
    'Roll disturbance (N m)',           {'env_bus.roll_dist'}
    });
G(end+1) = grp('navigation', {
    'Position: estimate vs truth (m)',  {'nav_bus.x_nav:i1-3', 'dyn_bus.pos_ned'}
    'Velocity: estimate vs truth (m/s)', {'nav_bus.x_nav:i4-6', 'dyn_bus.vel_ned'}
    'Filter diagnostics',               {'nav_bus.nav_diag'}
    'Shadow estimator',                 {'nav_bus.shadow'}
    });
G(end+1) = grp('fdir', {
    'Nav sensor mask',                  {'gnc_bus.nav_mask'}
    'Nav R scale',                      {'gnc_bus.nav_r_scale'}
    'Sensor valid flags',               {'nav_bus.valid'}
    });
end

function g = grp(name, displays)
sig = [displays{:, 2}];
b = cellfun(@(s) extractBefore(s, '.'), sig, 'uni', 0);
g = struct('name', name, 'buses', {unique(b(~strcmp(b, 'ref')), 'stable')}, 'displays', {displays});
end

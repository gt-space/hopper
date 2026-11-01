function [x_ref, a] = fdir_abort_reference(x, t, active, p, a)
% FDIR_ABORT_REFERENCE  R2 abort-to-land reference (GNC PDR slides 79, 84).
%
%   x       13-state [pos_ned; vel_ned; body rates; quaternion]
%   active  abort mode in effect this frame
%   p       .v_land (m/s sink at touchdown)  .v_max (m/s fastest abort
%           sink)  .tau (s, sink = height / tau above the floor)
%           .a_ramp (m/s^2, vertical ramp rate)  .dt (s, frame)
%   a       state (fdir_abort_state), returned updated
%
% Entry captures the state. The vertical velocity then ramps from the
% captured v_z toward a sink that shrinks with height (between v_land and
% v_max); the height reference integrates it. The horizontal command is
% zero velocity: the position reference follows the vehicle, so only
% velocity is regulated. Attitude goes to vertical, rates to zero. A
% fixed-gain regulator (fdir_command_select) flies it.
%
%#codegen

if ~active
    a.on = false;
    x_ref = x;
    return
end
if ~a.on
    a.on = true;
    a.t0 = t;
    a.vz = x(6);              % NED, down +
    a.z  = x(3);
end
alt = -a.z;
v_tgt = min(p.v_max, max(p.v_land, alt / p.tau));
step = p.a_ramp * p.dt;
a.vz = a.vz + min(max(v_tgt - a.vz, -step), step);
a.z = min(a.z + a.vz * p.dt, 0);

x_ref = zeros(13, 1);
x_ref(1:2) = x(1:2);                          % horizontal: zero-velocity hold
x_ref(3) = a.z;
x_ref(6) = a.vz;
x_ref(10:13) = [sqrt(2)/2; 0; sqrt(2)/2; 0];  % body x up
end

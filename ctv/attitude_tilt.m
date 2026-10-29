function tilt_deg = attitude_tilt(pitch_deg, yaw_deg)
%ATTITUDE_TILT Angle between the vehicle's body x axis and local up, in degrees.
%   Body x points up, so an upright vehicle logs pitch = 90 deg by construction
%   and the roll and yaw channels are undefined there. The direction cosine
%   matrix stays well defined, so tilt is the attitude metric to use instead of
%   max |pitch| or max |roll|.
%
%   Body x in NED is the first column of Rz(yaw)*Ry(pitch)*Rx(roll). Roll drops
%   out of that column, which is exactly why tilt survives the degeneracy.

th  = deg2rad(pitch_deg(:));
psi = deg2rad(yaw_deg(:));

bx = [cos(th) .* cos(psi), cos(th) .* sin(psi), -sin(th)];
up = [0 0 -1];                      % local up in NED, where z points down

tilt_deg = rad2deg(acos(max(-1, min(1, bx * up.'))));
end

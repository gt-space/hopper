function q = nav_quat_from_dcm(A)
%NAV_QUAT_FROM_DCM  Scalar-first quaternion of an inertial-to-body DCM.
%   q = nav_quat_from_dcm(A)   inverse of A = nav_dcm_from_quat(q)
%
%   A(q) = (w^2 - |v|^2) I + 2 v v' - 2 w [v x]: the attitude matrix the
%   nav filter uses (quaternion_to_rotation returns its transpose, body to
%   inertial) and the one the sim's 6DOF block uses (quat2dcm, NED to
%   body). So the sim and nav quaternions of their own DCMs share this
%   function. Shepperd's method (largest pivot), with w >= 0: the hopper
%   stays within ~30 deg of its pad attitude, a 90 deg rotation from NED,
%   so w never nears zero and the sign is continuous.
%
%#codegen

t = trace(A);
d = [t; A(1,1); A(2,2); A(3,3)];
[~, k] = max(d);
switch k
    case 1
        w = 0.5 * sqrt(1 + t);
        q = [w; (A(2,3) - A(3,2)) / (4*w); (A(3,1) - A(1,3)) / (4*w); (A(1,2) - A(2,1)) / (4*w)];
    case 2
        x = 0.5 * sqrt(1 + 2*A(1,1) - t);
        q = [(A(2,3) - A(3,2)) / (4*x); x; (A(1,2) + A(2,1)) / (4*x); (A(1,3) + A(3,1)) / (4*x)];
    case 3
        y = 0.5 * sqrt(1 + 2*A(2,2) - t);
        q = [(A(3,1) - A(1,3)) / (4*y); (A(1,2) + A(2,1)) / (4*y); y; (A(2,3) + A(3,2)) / (4*y)];
    otherwise
        z = 0.5 * sqrt(1 + 2*A(3,3) - t);
        q = [(A(1,2) - A(2,1)) / (4*z); (A(1,3) + A(3,1)) / (4*z); (A(2,3) + A(3,2)) / (4*z); z];
end
if q(1) < 0, q = -q; end
q = q / norm(q);
end

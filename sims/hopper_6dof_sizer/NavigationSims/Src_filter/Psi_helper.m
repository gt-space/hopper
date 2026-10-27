function Psi_out = Psi_helper(q)
% Maps an inertial-frame attitude error to a quaternion increment: [0;e] * q
q = q(:); w = q(1); x = q(2); y = q(3); z = q(4);
Psi_out = [ -x, -y, -z;
    w,  z, -y;
    -z,  w,  x;
    y, -x,  w ];
end
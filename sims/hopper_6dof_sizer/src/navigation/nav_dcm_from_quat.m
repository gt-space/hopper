function A = nav_dcm_from_quat(q)
%NAV_DCM_FROM_QUAT  Inertial-to-body DCM of a scalar-first quaternion.
%   A = nav_dcm_from_quat(q)   the transpose of the nav filter's
%                              quaternion_to_rotation (body to inertial)
%
%#codegen

A = quaternion_to_rotation(q)';
end

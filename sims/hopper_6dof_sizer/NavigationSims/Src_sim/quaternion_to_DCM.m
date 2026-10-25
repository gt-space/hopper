function rotationBodyToInertial = quaternion_to_DCM(quaternion)
% (Yuma's implementation | 10/24/26)
% Returns the rotation matrix from body to inertial frame associated with a quaternion
%
% Inputs:
% 
%   quaternion  - scalar-first quaternion (4x1)
% 
% Outputs:
%
%   rotationBodyToInertial - 3x3 rotation matrix transforming body-fixed vectors to inertial frame (3x3)

    assert(numel(quaternion) == 4, 'Input must have 4 elements'); % checks for quaternion size
    
    quaternion = quaternion(:); % ensures column vector format
    
    w = quaternion(1);
    x = quaternion(2);
    y = quaternion(3);
    z = quaternion(4);
    
    rotationBodyToInertial = ...
        [   1-2*(y^2+z^2),      2*(x*y+w*z),    2*(x*z-w*y);
              2*(x*y-w*z),    1-2*(x^2+z^2),    2*(y*z+w*x);
              2*(x*z+w*y),      2*(y*z-w*x),  1-2*(x^2+y^2)     ]';
end
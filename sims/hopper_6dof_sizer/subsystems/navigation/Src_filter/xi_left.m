function xi_matrix = xi_left(quaternion)
% (Yuma's implementation | 10/24/26)
% Returns the left quaternion multiplication / mapping matrix Xi_L(q)
%
% Inputs:
% 
%   quaternion  - Some scalar-first quaternion (4x1)
% 
% Outputs:
%
%   xi_matrix   - Left-multiplication xi matrix (4x3)
    
    assert(numel(quaternion) == 4, 'Input must have 4 elements'); % checks for quaternion size
    
    quaternion = quaternion(:); % ensures column vector format
    
    w = quaternion(1);
    x = quaternion(2);
    y = quaternion(3);
    z = quaternion(4);
    
    xi_matrix = ...
    [   -x,     -y,     -z;
         w,      z,     -y;
        -z,      w,      x;
         y,     -x,      w     ];
end
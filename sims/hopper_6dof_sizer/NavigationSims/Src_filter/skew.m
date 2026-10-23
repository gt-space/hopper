function S = skew(v)
% (Yuma's implementation | 10/23/26)
% Returns skew-symmetric matrix associated with a 3D vector
%
% Inputs:
% 
%   v   - Some generic vector (3x1)
% 
% Outputs:
%
%   S   - Skew-symmetric matrix satisfying S*u = cross(v,u) (3x3)

    assert(numel(v) == 3, 'Input must have 3 elements'); % checks for vector size

    v = v(:); % ensures column vector format

    S = [   0,  -v(3),  v(2);
         v(3),      0, -v(1);
        -v(2),   v(1),     0; ];
end


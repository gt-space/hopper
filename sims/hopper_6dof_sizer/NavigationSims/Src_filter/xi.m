function Xi = xi(q)
% (Yuma's implementation | 10/24/26)
% Returns the quaternion mapping matrix Xi associated with a 4x1 quaternion
%
% Inputs:
% 
%   q   - 4x1 quaternion (scalar first)
% 
% Outputs:
%
%   Xi  - 4x3 quaternion mapping matrix (4x3)

    assert(numel(q) == 4, 'Input must have 4 elements'); % checks for quaternion size

    q = q(:); % ensures column vector format

    w = q(1); x = q(2); y = q(3); z = q(4);

    Xi = [  -x, -y, -z;
             w, -z,  y;
             z,  w, -x;
            -y,  x,  w  ];
end
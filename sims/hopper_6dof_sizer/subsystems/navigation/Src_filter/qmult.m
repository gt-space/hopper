%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qmult.m
% 
% Author: David Meisinger
% Date: 9/19/2026
% 
% Takes two quaternions and returns the multiplication of them. Assumes
% passive rotation convention. Can accommodate lists of quaternions.
% 
% Inputs:
% q1 - 4xn quaternion (scalar first)
% q2 - 4xn quaternion (scalar first)
% 
% Outputs:
% q - 4xn quaternion, equals q1 multiplied by q2
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q = qmult(q1, q2)

    % Validate sizes
    assert(size(q1,1)==4, 'q1 must be 4×n.');
    assert(size(q2,1)==4, 'q2 must be 4xn.');
    assert(size(q1,2)==size(q2,2), 'Quaternion lists must be of same length.');
    
    % Separate the scalar and vector parts
    v1 = q1(2:4,:);
    s1 = q1(1,:);
    v2 = q2(2:4,:);
    s2 = q2(1,:);

    % Perform the multiplication
    q = [s1.*s2 - dot(v1, v2, 1); (v1.*s2 + v2.*s1 + cross(v1, v2, 1))];
    
end
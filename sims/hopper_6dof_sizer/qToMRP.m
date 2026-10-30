function [p] = qToMRP(q)
%input scalar first quaternion

qv = q(2:4);
qs = q(1);

p = qv/(1+qs);

end
function [Atilde,Btilde] = extendedDynamics(A,B)

Atilde = zeros(15,15);
Btilde = zeros(15,4);

Atilde(1:12,1:12) = A;
Atilde(13:15,:) = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                   0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
                   0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
Btilde(1:12,1:4) =B;
Btilde(13:15,:) = [0,0,0,0;
                   0,0,0,0;
                   0,0,0,0];

end
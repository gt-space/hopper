function [Aa,Ba,Cpos,Dpos] = extendedDynamics(A,B)

    nx = size(A,1);
    nu = size(B,2);

    Cpos = [eye(3), zeros(3,nx-3)];
    Dpos = zeros(3,nu);

    % xi_dot = r - y
    Aa = [ A,     zeros(nx,3);
          -Cpos,  zeros(3,3)];

    Ba = [ B;
          -Dpos];
end
function Vdot = Riccati_V(~, V, Ai, Bi, Si, R_riccati, xi_ref)
    % Vdot = A*x_ref - A'*V + S*B*R^{-1}*B'*V
    R_inv = inv(R_riccati);
    Vdot = Ai*xi_ref - Ai'*V + Si*Bi*R_inv*Bi'*V;
end
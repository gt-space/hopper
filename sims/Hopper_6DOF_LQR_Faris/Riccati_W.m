function Wdot = Riccati_W(~, ~, Bi, Q_riccati, R_riccati, xi_ref, Vi)
    % Wdot = -0.5*x_ref'*Q*x_ref + 0.5*V'*B*R^{-1}*B'*V
    R_inv = inv(R_riccati);
    Wdot = -0.5 * xi_ref' * Q_riccati * xi_ref + 0.5 * Vi' * Bi * R_inv * Bi' * Vi;
end
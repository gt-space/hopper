function Sdot = Riccati(~, S_vec, Ai, Bi, Q_riccati, R_riccati)
    % Reshape vector to matrix
    S_mat = reshape(S_vec, [13, 13]);
    
    % Riccati equation: dS/dt = -A'*S - S*A - Q + S*B*R^{-1}*B'*S
    R_inv = inv(R_riccati);
    Sdot_mat = -Ai'*S_mat - S_mat*Ai - Q_riccati + S_mat*Bi*R_inv*Bi'*S_mat;
    
    % Return as vector
    Sdot = Sdot_mat(:);
end
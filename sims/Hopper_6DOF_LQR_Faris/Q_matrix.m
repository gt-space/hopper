function   Qdot = Q_matrix(~, R_current, ~, Bi, ~, R_riccati )
    Qdot_mat = R_current' * Bi * inv(R_riccati) * Bi' * R_current;    
    Qdot = Qdot_mat(:);                % return as vector for ode45
end

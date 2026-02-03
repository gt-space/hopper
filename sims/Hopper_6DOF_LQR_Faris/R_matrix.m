function   Rdot = R_matrix(~, S_current, R, Ai, Bi, ~, R_riccati )

    Rmat = reshape(R, [13 13]);        % convert state vector to 13x13 matrix
    Rdot_mat = - ( Ai' - S_current * Bi * (R_riccati \ Bi') ) * Rmat;    
    Rdot = Rdot_mat(:);                % return as vector for ode45
end


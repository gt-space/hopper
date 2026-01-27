function   Sdot = Riccati(~, S, Ai, Bi, Q_riccati, R_riccati)

  Smat = reshape(S, [13 13]);        % convert state vector to 13x13 matrix
  temp = R_riccati \ (Bi' * Smat);   
  Sdot_mat = -Smat * Ai - Ai' * Smat - Q_riccati + Smat * Bi * temp;
  Sdot = Sdot_mat(:);                % return as vector for ode45
end

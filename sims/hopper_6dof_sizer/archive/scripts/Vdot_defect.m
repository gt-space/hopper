function dV = Vdot_defect(t, V, tS, S_hist, A, B, Q, R, d, tgrid)

nx = size(d,2);

% interpolate A(t), B(t), d(t), S(t)
Ai = interp3D(A, tgrid, t);
Bi = interp3D(B, tgrid, t);
di = interpVec(d, tgrid, t);          % nx x 1
Si = interp3D(S_hist, tS, t);

% defect-based V dynamics
dV = -(Ai' - Si*Bi*(R\Bi')) * V - Si * di;

dV = dV(:);

end
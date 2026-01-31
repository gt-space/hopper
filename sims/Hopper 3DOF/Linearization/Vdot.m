function dV = Vdot(t,V,tS,Ssol,Adyn,Bdyn,Q,R,eta,tgrid)
    Ai = lookupMat(tgrid, Adyn, t);
    Bi = lookupMat(tgrid, Bdyn, t);

    % interpolate S using its OWN time vector tS
    Si = lookupMat(tS, Ssol, t);

    eta_t = interp1(tgrid, eta, t)';   % eta is N×6

    dV = -(-Q*eta_t + (Ai' - Si*Bi*(R\Bi'))*V);
end
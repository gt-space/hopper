function dS = Sdot(t,Svec,Adyn,Bdyn,Q,R,tgrid)
    nx = 6;
    S = reshape(Svec,nx,nx);
    
    A = lookupMat(tgrid, Adyn, t);
    B = lookupMat(tgrid, Bdyn, t);

    
    dSmat = -(Q + S*A + A'*S - S*B*(R\B')*S);
    dS = dSmat(:);
end
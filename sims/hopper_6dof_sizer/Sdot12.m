function dS = Sdot12(t,Svec,Adyn,Bdyn,Q,R,tgrid)
    nx = 12;
    S = reshape(Svec,nx,nx);
    
    A = lookupMat(tgrid, Adyn, t);
    B = lookupMat(tgrid, Bdyn, t);

    
    dSmat = -(Q + S*A + A'*S - S*B*(R\B')*S);
    %dSmat = -(S*A+A'*S+Q-S*B*R^-1*B'*S);
    dS = dSmat(:);
end
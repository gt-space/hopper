function vi = interpVec(vhist, tgrid, t)
% vhist is N x n
vi = interp1(tgrid, vhist, t, 'linear', 'extrap').';
end
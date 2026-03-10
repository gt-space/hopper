function Mi = interp3D(Mhist, tgrid, t)

[n1,n2,N] = size(Mhist);

Mflat = reshape(Mhist, n1*n2, N).';             % N x (n1*n2)
Mflat_i = interp1(tgrid, Mflat, t, 'linear', 'extrap');   % 1 x (n1*n2)
Mi = reshape(Mflat_i.', n1, n2);

end
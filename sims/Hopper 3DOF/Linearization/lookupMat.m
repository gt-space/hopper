function Ak = lookupMat(tgrid, Ahist, t)
    % Ahist: n×m×N
    [~,k] = min(abs(tgrid - t));   % nearest index
    Ak = Ahist(:,:,k);
end
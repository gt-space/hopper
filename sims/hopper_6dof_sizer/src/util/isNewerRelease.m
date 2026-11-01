function tf = isNewerRelease(rel, than)
%ISNEWERRELEASE  True if MATLAB release rel is later than release than.
%   isNewerRelease('R2026b', 'R2026a')   % true
%   isNewerRelease('R2022b', teamRelease())

tf = releaseOrder(rel) > releaseOrder(than);
end

function n = releaseOrder(rel)
% 'R2025b' -> 4051
n = 2 * str2double(rel(2:5)) + (lower(rel(6)) == 'b');
end

function nw = mcWorkerCount(nRuns, perWorkerGB, reserveGB)
%MCWORKERCOUNT  Parallel workers that fit in the RAM available right now.
%   nw = mcWorkerCount(nRuns) is at most nRuns, the number of physical
%   cores, and (free RAM - reserve) / per-worker memory, and at least 1.
%   Each worker holds MATLAB, Simulink and the compiled 6DOF model, about
%   3.3 GB committed (measured on R2026a, 2026-10-24). Asking for more than
%   fits makes the OS page to disk; on Windows the page file can grow by
%   tens of GB.

if nargin < 2, perWorkerGB = 3.5; end
if nargin < 3, reserveGB = 4; end

freeGB = availableMemoryGB();
cores = feature('numcores');
nw = max(1, min([nRuns, cores, floor((freeGB - reserveGB) / perWorkerGB)]));
end

function gb = availableMemoryGB()
if ispc
    [~, sys] = memory();
    gb = sys.PhysicalMemory.Available / 2^30;
elseif isfile('/proc/meminfo')
    txt = fileread('/proc/meminfo');
    tok = regexp(txt, 'MemAvailable:\s*(\d+)', 'tokens', 'once');
    gb = str2double(tok{1}) / 2^20;   % kB -> GB
else
    gb = 8;   % unknown platform: assume a small machine
end
end

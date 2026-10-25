function outDir = export_gains(outDir)
%EXPORT_GAINS  Write the LQR gain schedule to CSV for the Rust controller.
%   export_gains() writes tgrid, xref (Target_Trajectory), k1flat, k2grid and
%   unom from the base workspace to fsw_bridge/build/gains/*.csv, one row per
%   grid time, full double precision. These are the same tables the
%   Simulink Controller block reads (Constant4..8).

here = fileparts(mfilename('fullpath'));
if nargin < 1, outDir = fullfile(here, 'build', 'gains'); end
if ~isfolder(outDir), mkdir(outDir); end

if ~evalin('base', 'exist(''tgrid'', ''var'') && exist(''K1flat'', ''var'')')
    cd(fileparts(here));
    evalin('base', 'sim_setup_cached');
end

tables = {'tgrid', 'tgrid'; 'xref', 'Target_Trajectory'; 'k1flat', 'K1flat'; ...
          'k2grid', 'K2grid'; 'unom', 'unom'};
for i = 1:size(tables, 1)
    A = evalin('base', tables{i, 2});
    if isvector(A), A = A(:); end
    writematrix(A, fullfile(outDir, [tables{i, 1} '.csv']), 'Delimiter', ',');
    % writematrix keeps 15 significant digits; rewrite at full precision
    fid = fopen(fullfile(outDir, [tables{i, 1} '.csv']), 'w');
    fmt = [strjoin(repmat({'%.17g'}, 1, size(A, 2)), ',') '\n'];
    fprintf(fid, fmt, A.');
    fclose(fid);
    fprintf('%-7s %-18s %d x %d\n', tables{i, 1}, tables{i, 2}, size(A, 1), size(A, 2));
end
end

function reports = modelDiffReport(baseRef, outDir)
%MODELDIFFREPORT  HTML Simulink comparison of every .slx changed since baseRef.
%   reports = modelDiffReport('origin/main')
%   reports = modelDiffReport(baseRef, outDir)   default outDir outputs/model-diffs
%
%   For each model under this project that differs between baseRef and the
%   working tree, the baseRef version is extracted with git show and
%   compared with visdiff; the report is published as HTML. CI uploads the
%   folder as an artifact on pull requests so reviewers can see what changed
%   inside binary models.
%
%   Needs only tools/ci and src/util on the path, not the open project: it
%   compares two files. Every step is also appended, with a timestamp, to
%   progress.log in outDir and flushed at once, so if the run stalls or is
%   killed the uploaded log shows how far it got.

if nargin < 2, outDir = hopperRoot('outputs', 'model-diffs'); end
if ~isfolder(outDir), mkdir(outDir); end
logFile = fullfile(outDir, 'progress.log');
note = @(fmt, varargin) progress(logFile, fmt, varargin{:});
note('start, base %s, MATLAB %s', baseRef, version);

root = hopperRoot();
[st, changed] = system(sprintf('git -C "%s" diff --name-only --relative "%s" -- "*.slx"', root, baseRef));
assert(st == 0, 'git diff failed: %s', changed);
changed = splitlines(strtrim(changed));
changed = changed(~cellfun(@isempty, changed));
note('%d changed model(s)', numel(changed));

reports = {};
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
for k = 1:numel(changed)
    rel = changed{k};
    newFile = fullfile(root, rel);
    [~, name] = fileparts(rel);
    if ~isfile(newFile)
        note('%s: deleted', rel); continue
    end
    oldFile = fullfile(tmp, [name '_base.slx']);
    st = system(sprintf('git -C "%s" show "%s:./%s" > "%s"', root, baseRef, strrep(rel, '\', '/'), oldFile));
    if st ~= 0
        note('%s: new file', rel); continue
    end
    try
        note('%s: comparing', rel);
        comp = visdiff(oldFile, newFile);
        note('%s: publishing', rel);
        html = publish(comp, Format='html', OutputFolder=outDir, ...
            Name=matlab.lang.makeValidName(erase(rel, '.slx')));
        reports{end+1} = html; %#ok<AGROW>
        note('%s: %s', rel, html);
    catch e
        note('%s: comparison failed (%s)', rel, e.message);
    end
    bdclose('all');
end
note('%d model diff report(s) in %s', numel(reports), outDir);
end

function progress(logFile, fmt, varargin)
% to the console and, opened and closed each time, to the log file
msg = sprintf(['%s  ' fmt], char(datetime('now', Format='HH:mm:ss')), varargin{:});
fprintf('%s\n', msg);
fid = fopen(logFile, 'a');
if fid > 0
    fprintf(fid, '%s\n', msg);
    fclose(fid);
end
end

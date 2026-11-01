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

if nargin < 2, outDir = hopperRoot('outputs', 'model-diffs'); end
if ~isfolder(outDir), mkdir(outDir); end
root = hopperRoot();
[st, changed] = system(sprintf('git -C "%s" diff --name-only --relative "%s" -- "*.slx"', root, baseRef));
assert(st == 0, 'git diff failed: %s', changed);
changed = splitlines(strtrim(changed));
changed = changed(~cellfun(@isempty, changed));

reports = {};
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
for k = 1:numel(changed)
    rel = changed{k};
    newFile = fullfile(root, rel);
    [~, name] = fileparts(rel);
    if ~isfile(newFile)
        fprintf('%s: deleted\n', rel); continue
    end
    oldFile = fullfile(tmp, [name '_base.slx']);
    st = system(sprintf('git -C "%s" show "%s:./%s" > "%s"', root, baseRef, strrep(rel, '\', '/'), oldFile));
    if st ~= 0
        fprintf('%s: new file\n', rel); continue
    end
    try
        comp = visdiff(oldFile, newFile);
        html = publish(comp, Format='html', OutputFolder=outDir, ...
            Name=matlab.lang.makeValidName(erase(rel, '.slx')));
        reports{end+1} = html; %#ok<AGROW>
        fprintf('%s: %s\n', rel, html);
    catch e
        fprintf('%s: comparison failed (%s)\n', rel, e.message);
    end
    bdclose('all');
end
fprintf('%d model diff report(s) in %s\n', numel(reports), outDir);
end

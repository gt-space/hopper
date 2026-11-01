function syncProjectFiles()
%SYNCPROJECTFILES  Make the HopperSim project file list match the folders.
%   Run after moving or adding files outside the MATLAB Project UI (for
%   example with git mv). Files that no longer exist are removed from the
%   project; new files under the project's folders are added. Search-path
%   folders are not changed; add a new folder with addPath(currentProject, f).

proj = currentProject();
root = proj.RootFolder;

listed = [proj.Files.Path];
gone = listed(~isfile(listed) & ~isfolder(listed));
for k = 1:numel(gone)
    removeFile(proj, gone(k));
    fprintf('removed  %s\n', erase(gone(k), root + filesep));
end

ignore = ["outputs", "work", "resources", "slprj", ".git"];
f = dir(fullfile(root, '**', '*'));
f = f(~[f.isdir]);
paths = string(fullfile({f.folder}, {f.name}));
rel = erase(paths, root + filesep);
top = extractBefore(rel + filesep, filesep);
keep = ~ismember(top, ignore) & ~endsWith(rel, [".asv", ".slxc", ".autosave", ".prj"]) & ...
    ~startsWith(rel, ".");
listed = [proj.Files.Path];
for p = paths(keep & ~ismember(paths, listed))
    try
        addFile(proj, p);
        fprintf('added    %s\n', erase(p, root + filesep));
    catch e
        fprintf('skipped  %s (%s)\n', erase(p, root + filesep), e.message);
    end
end
end

function rel = teamRelease()
%TEAMRELEASE  The MATLAB release the team saves models in, e.g. 'R2026a'.
%   Read from the TEAM_RELEASE file at the project root, which is also what
%   CI installs and what tools/ci/guards.sh checks. To move the team to a
%   new release, edit that one line (see docs/CONTRIBUTING.md, "Moving to a
%   new MATLAB release").

rel = strtrim(fileread(hopperRoot('TEAM_RELEASE')));
assert(~isempty(regexp(rel, '^R\d{4}[ab]$', 'once')), ...
    'TEAM_RELEASE must hold a release name like R2026a, not "%s".', rel);
end

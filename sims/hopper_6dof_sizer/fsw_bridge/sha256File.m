function h = sha256File(f)
%SHA256FILE  Lower-case hex SHA-256 of a file.
%   Uses the operating system's tool (certutil on Windows, sha256sum or
%   shasum elsewhere), so it works in matlab -batch without a Java runtime.

f = char(f);
assert(isfile(f), 'sha256File:missing', '%s does not exist', f);
if ispc
    [s, r] = system(sprintf('certutil -hashfile "%s" SHA256', f));
    tok = regexp(r, '^\s*([0-9a-fA-F ]{64,95})\s*$', 'tokens', 'lineanchors', 'once');
else
    [s, r] = system(sprintf('sha256sum "%s" 2>/dev/null || shasum -a 256 "%s"', f, f));
    tok = regexp(r, '^([0-9a-fA-F]{64})', 'tokens', 'once');
end
assert(s == 0 && ~isempty(tok), 'sha256File:tool', 'could not hash %s: %s', f, r);
h = lower(strrep(tok{1}, ' ', ''));
assert(numel(h) == 64, 'sha256File:format', 'unexpected hash output for %s: %s', f, r);
end

"""Fail if GitHub can't resolve an owner in .github/CODEOWNERS.

    python tools/ci/codeowners_check.py <owner/repo> <ref> [token]

GitHub silently ignores CODEOWNERS lines whose team or user doesn't exist or
lacks write access, so a renamed or deleted team would just stop requesting
reviews. This asks GitHub's codeowners/errors API for the given commit and
turns each error into an ::error annotation.
"""
import json
import sys
import urllib.error
import urllib.request


def main(repo, ref, token=None):
    url = f"https://api.github.com/repos/{repo}/codeowners/errors?ref={ref}"
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json"})
    if token:
        req.add_header("Authorization", f"Bearer {token}")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            errors = json.load(resp).get("errors", [])
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print("CODEOWNERS: no CODEOWNERS file at this commit; nothing to check")
            return
        raise
    for e in errors:
        print(f"::error file={e.get('path', '.github/CODEOWNERS')},line={e.get('line', 1)}::"
              f"{e.get('kind', 'CODEOWNERS error')}: {e.get('suggestion') or e.get('message', '').splitlines()[0]}")
    if errors:
        print(f"CODEOWNERS: {len(errors)} line(s) GitHub can't resolve; see the annotations. "
              "A team needs Write access to the repo to own code.")
        sys.exit(1)
    print("CODEOWNERS: every owner resolves")


if __name__ == "__main__":
    main(*sys.argv[1:])

"""Turn buildtool's JUnit XML into GitHub annotations and a job-summary table.

    python tools/ci/junit_report.py outputs/test-results/*.xml >> "$GITHUB_STEP_SUMMARY"

Each failed test becomes an ::error annotation, which shows on the run page and
in a PR's checks without signing in. The summary gets one row per failed test.
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def clean(text, limit):
    text = " ".join((text or "").split())
    return text if len(text) <= limit else text[: limit - 3] + "..."


def main(paths):
    total = failed = 0
    rows = []
    for path in paths:
        p = Path(path)
        if not p.is_file():
            continue
        for case in ET.parse(p).getroot().iter("testcase"):
            total += 1
            bad = case.find("failure")
            if bad is None:
                bad = case.find("error")
            if bad is None:
                continue
            failed += 1
            name = f"{case.get('classname', '')}/{case.get('name', '')}"
            detail = bad.get("message") or bad.text or ""
            # stderr: annotation; stdout: summary markdown
            print(f"::error title={clean(name, 120)}::{clean(detail, 1500)}", file=sys.stderr)
            rows.append((p.stem, name, clean(detail, 400).replace("|", "\\|")))

    print("## Test results\n")
    if total == 0:
        print("No test results were written; the build stopped before running tests. See the job log.")
        return
    print(f"**{total - failed}/{total}** passed.\n")
    if rows:
        print("| Task | Test | Failure |\n|---|---|---|")
        for task, name, detail in rows:
            print(f"| {task} | `{name}` | {detail} |")


if __name__ == "__main__":
    main(sys.argv[1:])

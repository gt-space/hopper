"""Render the Monte Carlo summary CSV (runMonteCarlo) as a GitHub job-summary table.

    python tools/ci/mc_summary.py outputs/montecarlo/mc_summary.csv >> "$GITHUB_STEP_SUMMARY"
"""
import csv
import sys
from pathlib import Path

BOOL = "bool"
COLUMNS = [
    ("run", "Run", None),
    ("pass", "CI limits", BOOL),
    ("max_alt", "Max alt (m)", 2),
    ("landing_vel", "Impact (m/s)", 3),
    ("landing_offset", "Offset (m)", 2),
    ("max_tilt", "Max tilt (deg)", 2),
    ("flight_time", "Flight (s)", 2),
    ("meets_requirements", "Mission req.", BOOL),
    ("missed_requirements", "Missed requirements", None),
]


def cell(value, digits):
    if digits == BOOL:
        return "pass" if value in ("true", "1") else "**fail**"
    if digits is not None:
        try:
            return f"{float(value):.{digits}f}"
        except ValueError:
            return value
    return value


def main(path):
    p = Path(path)
    print("## Monte Carlo (10 seeded runs)\n")
    if not p.is_file():
        print("No summary was written; the Monte Carlo did not finish. See the job log.")
        return
    rows = list(csv.DictReader(p.open(newline="")))
    passed = sum(r["pass"] in ("true", "1") for r in rows)
    met = sum(r["meets_requirements"] in ("true", "1") for r in rows)
    print(f"**{passed}/{len(rows)}** within CI regression limits (`tests/montecarlo/mcLimits.m`), "
          f"**{met}/{len(rows)}** meet the mission requirements (`inputs/mission_inputs.m`).\n")
    print("| " + " | ".join(h for _, h, _ in COLUMNS) + " |")
    print("|" + "---|" * len(COLUMNS))
    for r in rows:
        print("| " + " | ".join(cell(r.get(k, ""), d) for k, _, d in COLUMNS) + " |")


if __name__ == "__main__":
    main(sys.argv[1])

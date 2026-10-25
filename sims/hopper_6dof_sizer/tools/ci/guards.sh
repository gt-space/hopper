#!/usr/bin/env bash
# Repo guards for the HopperSim project that need no MATLAB (CI "guards" job).
# The same checks run inside MATLAB as `buildtool check`.
#   - no Simulink build caches, autosaves or .slx backups are committed
#   - no live model (outside archive/) is saved newer than the team release
# Run from sims/hopper_6dof_sizer.
set -euo pipefail

TEAM_RELEASE="R2026a"
fail=0

rank() {  # R2025b -> 4051
  local r="$1"; local y="${r:1:4}"; local b=0
  [[ "${r:5:1}" == "b" ]] && b=1
  echo $(( 2 * 10#$y + b ))
}

junk=$(git ls-files | grep -E '(^|/)slprj/|\.slxc$|\.asv$|\.autosave$|\.slx\.r20[0-9]{2}[ab]$' || true)
if [[ -n "$junk" ]]; then
  echo "::error::Build caches / autosaves are committed (they are gitignored; remove them):"
  echo "$junk" | sed 's/^/  /'
  fail=1
fi

team=$(rank "$TEAM_RELEASE")
newer=""
while IFS= read -r f; do
  rel=$(unzip -p "$f" metadata/coreProperties.xml 2>/dev/null | grep -o '<cp:version>R[0-9]\{4\}[ab]</cp:version>' | sed 's/<[^>]*>//g' || true)
  if [[ -z "$rel" ]]; then
    echo "::warning file=sims/hopper_6dof_sizer/$f::could not read the saved release"
    continue
  fi
  if (( $(rank "$rel") > team )); then
    newer+="  $f ($rel)"$'\n'
    echo "::error file=sims/hopper_6dof_sizer/$f::saved in $rel; export it to $TEAM_RELEASE (Simulink.exportToVersion)"
  fi
done < <(git ls-files '*.slx' | grep -v '^archive/')
if [[ -n "$newer" ]]; then
  echo "Models saved newer than $TEAM_RELEASE:"; printf '%s' "$newer"
  fail=1
fi

n=$(git ls-files '*.slx' | grep -vc '^archive/' || true)
if (( fail == 0 )); then
  echo "guards: OK ($n live models, none newer than $TEAM_RELEASE, no committed caches)"
fi
exit $fail

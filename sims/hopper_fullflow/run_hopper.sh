#!/usr/bin/env bash
# Run from Git Bash without relying on an editor-generated Windows drive path.
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd -- "$project_dir"

for environment in .venv .venv-validation venv; do
    interpreter="./$environment/Scripts/python.exe"
    if [[ -f "$interpreter" ]] && "$interpreter" -c 'import fullflow, fullplot, thermoprop, numpy, scipy, h5py, matplotlib' >/dev/null 2>&1; then
        exec "$interpreter" Hopper.py "$@"
    fi
done

printf '%s\n' 'No working project Python environment with the required dependencies was found.' >&2
printf '%s\n' 'Create a Windows Python venv in .venv and install requirements.txt.' >&2
exit 1

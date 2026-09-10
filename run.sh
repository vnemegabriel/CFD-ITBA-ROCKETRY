#!/bin/bash
# ./run.sh [sweep.txt] [runs-dir]        default: sweep.txt, ./runs
#
# Builds and runs every case in the sweep file, one after the other.
# A case whose directory already exists is skipped, so re-running resumes.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
SWEEP=${1:-$ROOT/sweep.txt}
RUNS=${2:-$ROOT/runs}

[ -f "$SWEEP" ] || { echo "no sweep file: $SWEEP"; exit 1; }
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
mkdir -p "$RUNS"

while read -r NAME OPTS; do
    case "$NAME" in ''|'#'*) continue ;; esac
    DIR="$RUNS/$NAME"

    if [ -d "$DIR" ]; then
        echo "== $NAME  already exists, skipping"
        continue
    fi

    echo "== $NAME  $OPTS"
    if ! "$ROOT/newCase.sh" "$DIR" $OPTS; then
        echo "!! $NAME failed to set up, stopping"
        exit 1
    fi

    ( cd "$DIR" && ./Allrun ) || echo "!! $NAME failed to run, continuing"
done < "$SWEEP"

echo
echo "done.  results under $RUNS"

#!/bin/bash
# ./runCase.sh [sweep.txt] [runs-dir] [meshes-dir]
#                                    default: sweep.txt, ~/runs, ~/meshes
#
# Makes and runs every case in the sweep file, one after the other.  Each line
# names its mesh with --mesh; a bare name lives in meshes-dir.  The meshes must
# exist: build them first with runMesh.sh.  A case whose directory already
# exists is skipped, so re-running resumes.
#
# The script reads the sweep file while it runs: do not edit it in place until
# it ends.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
SWEEP=${1:-$ROOT/sweep.txt}
RUNS=${2:-$HOME/runs}
MESHES=${3:-$HOME/meshes}

[ -f "$SWEEP" ] || { echo "no sweep file: $SWEEP"; exit 1; }
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
mkdir -p "$RUNS"

while read -r NAME REST; do
    case "$NAME" in ''|'#'*) continue ;; esac
    DIR="$RUNS/$NAME"

    if [ -d "$DIR" ]; then
        echo "== $NAME  already exists, skipping"
        continue
    fi

    MESH=""; OPTS=()
    set -- $REST
    while [ $# -gt 0 ]; do
        case "$1" in
            --mesh) MESH=$2; shift 2 ;;
            *)      OPTS+=("$1"); shift ;;
        esac
    done
    [ -n "$MESH" ] || { echo "!! $NAME has no --mesh, stopping"; exit 1; }
    MESH=${MESH/#\~/$HOME}
    case "$MESH" in */*) ;; *) MESH=$MESHES/$MESH ;; esac
    [ -d "$MESH/constant/polyMesh" ] || {
        echo "!! $NAME: no mesh in $MESH -- build it with runMesh.sh, stopping"; exit 1; }

    echo "== $NAME  --mesh $MESH ${OPTS[*]}"
    if ! "$ROOT/newCase.sh" "$DIR" --mesh "$MESH" "${OPTS[@]}" < /dev/null; then
        echo "!! $NAME failed to set up, stopping"
        exit 1
    fi

    ( cd "$DIR" && ./Allrun < /dev/null ) || echo "!! $NAME failed to run, continuing"
done < "$SWEEP"

echo
echo "done.  results under $RUNS"

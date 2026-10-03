#!/bin/bash
# ./run.sh [sweep.txt] [runs-dir] [meshes-dir]    default: sweep.txt, ~/runs, ~/meshes
#
# Builds and runs every case in the sweep file, one after the other.  Each line
# names its mesh with --mesh; a bare name lives in meshes-dir.  A mesh that does
# not exist yet is built first by mesh/Allmesh with the line's --maxCellSize and
# --finLevel, and reused by later lines.  A case whose directory already exists
# is skipped, so re-running resumes.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
SWEEP=${1:-$ROOT/sweep.txt}
RUNS=${2:-$HOME/runs}
MESHES=${3:-$HOME/meshes}

[ -f "$SWEEP" ] || { echo "no sweep file: $SWEEP"; exit 1; }
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
mkdir -p "$RUNS"

meshValue() { foamDictionary -entry "$2" -value "$1/system/meshDict" 2>/dev/null; }

while read -r NAME REST; do
    case "$NAME" in ''|'#'*) continue ;; esac
    DIR="$RUNS/$NAME"

    if [ -d "$DIR" ]; then
        echo "== $NAME  already exists, skipping"
        continue
    fi

    MESH=""; H=""; FIN=""; OPTS=()
    set -- $REST
    while [ $# -gt 0 ]; do
        case "$1" in
            --mesh)        MESH=$2; shift 2 ;;
            --maxCellSize) H=$2;    shift 2 ;;
            --finLevel)    FIN=$2;  shift 2 ;;
            *)             OPTS+=("$1"); shift ;;
        esac
    done
    [ -n "$MESH" ] || { echo "!! $NAME has no --mesh, stopping"; exit 1; }
    MESH=${MESH/#\~/$HOME}
    case "$MESH" in */*) ;; *) MESH=$MESHES/$MESH ;; esac

    if [ -d "$MESH/constant/polyMesh" ]; then
        if { [ -n "$H" ] && [ "$(meshValue "$MESH" maxCellSize)" != "$H" ]; } ||
           { [ -n "$FIN" ] && [ "$(meshValue "$MESH" localRefinement/fins/additionalRefinementLevels)" != "$FIN" ]; }; then
            echo "!! $NAME asks for a different $MESH than the one built, stopping"
            exit 1
        fi
    elif [ -e "$MESH" ]; then
        echo "!! $MESH exists but holds no mesh: delete it, stopping"
        exit 1
    else
        echo "== mesh $MESH"
        if ! "$ROOT/mesh/Allmesh" "$MESH" ${H:+--maxCellSize "$H"} ${FIN:+--finLevel "$FIN"} < /dev/null; then
            echo "!! $MESH failed to build, stopping"
            exit 1
        fi
    fi

    echo "== $NAME  --mesh $MESH ${OPTS[*]}"
    if ! "$ROOT/newCase.sh" "$DIR" --mesh "$MESH" "${OPTS[@]}" < /dev/null; then
        echo "!! $NAME failed to set up, stopping"
        exit 1
    fi

    ( cd "$DIR" && ./Allrun < /dev/null ) || echo "!! $NAME failed to run, continuing"
done < "$SWEEP"

echo
echo "done.  results under $RUNS"

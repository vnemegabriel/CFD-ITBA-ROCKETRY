#!/bin/bash
# ./run.sh [--meshOnly] [sweep.txt] [runs-dir] [meshes-dir]
#                                    default: sweep.txt, ~/runs, ~/meshes
#
# Builds and runs every case in the sweep file, one after the other.  Each line
# names its mesh with --mesh; a bare name lives in meshes-dir.  A mesh that does
# not exist yet is built first by mesh/Allmesh with the line's --maxCellSize,
# --finLevel, --nLayers and --thicknessRatio, and reused by later lines.  A case
# whose directory already exists is skipped, so re-running resumes.
# --meshOnly builds the meshes and stops there.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
MESHONLY=0
[ "${1:-}" = --meshOnly ] && { MESHONLY=1; shift; }
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

    if [ $MESHONLY = 0 ] && [ -d "$DIR" ]; then
        echo "== $NAME  already exists, skipping"
        continue
    fi

    MESH=""; H=""; FIN=""; NL=""; TR=""; OPTS=()
    set -- $REST
    while [ $# -gt 0 ]; do
        case "$1" in
            --mesh)        MESH=$2; shift 2 ;;
            --maxCellSize) H=$2;    shift 2 ;;
            --finLevel)    FIN=$2;  shift 2 ;;
            --nLayers)     NL=$2;   shift 2 ;;
            --thicknessRatio) TR=$2; shift 2 ;;
            *)             OPTS+=("$1"); shift ;;
        esac
    done
    [ -n "$MESH" ] || { echo "!! $NAME has no --mesh, stopping"; exit 1; }
    MESH=${MESH/#\~/$HOME}
    case "$MESH" in */*) ;; *) MESH=$MESHES/$MESH ;; esac

    if [ -d "$MESH/constant/polyMesh" ]; then
        if { [ -n "$H" ] && [ "$(meshValue "$MESH" maxCellSize)" != "$H" ]; } ||
           { [ -n "$FIN" ] && [ "$(meshValue "$MESH" localRefinement/fins/additionalRefinementLevels)" != "$FIN" ]; } ||
           { [ -n "$NL" ] && [ "$(meshValue "$MESH" boundaryLayers/patchBoundaryLayers/walls/nLayers)" != "$NL" ]; } ||
           { [ -n "$TR" ] && [ "$(meshValue "$MESH" boundaryLayers/patchBoundaryLayers/walls/thicknessRatio)" != "$TR" ]; }; then
            echo "!! $NAME asks for a different $MESH than the one built, stopping"
            exit 1
        fi
    elif [ -e "$MESH" ]; then
        echo "!! $MESH exists but holds no mesh: delete it, stopping"
        exit 1
    else
        echo "== mesh $MESH"
        if ! "$ROOT/mesh/Allmesh" "$MESH" ${H:+--maxCellSize "$H"} ${FIN:+--finLevel "$FIN"} \
                ${NL:+--nLayers "$NL"} ${TR:+--thicknessRatio "$TR"} < /dev/null; then
            echo "!! $MESH failed to build, stopping"
            exit 1
        fi
    fi

    [ $MESHONLY = 1 ] && continue

    echo "== $NAME  --mesh $MESH ${OPTS[*]}"
    if ! "$ROOT/newCase.sh" "$DIR" --mesh "$MESH" "${OPTS[@]}" < /dev/null; then
        echo "!! $NAME failed to set up, stopping"
        exit 1
    fi

    ( cd "$DIR" && ./Allrun < /dev/null ) || echo "!! $NAME failed to run, continuing"
done < "$SWEEP"

echo
if [ $MESHONLY = 1 ]; then echo "done.  meshes under $MESHES"; else echo "done.  results under $RUNS"; fi

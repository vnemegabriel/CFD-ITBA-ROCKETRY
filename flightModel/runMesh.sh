#!/bin/bash
# ./runMesh.sh [meshes.txt] [meshes-dir]
#                                    default: meshes.txt, ~/meshes
#
# Builds every mesh in the file, one after the other, with mesh/Allmesh.  Each
# line is a mesh name and the options of Allmesh.  A mesh that exists with the
# same options is kept; one built with different options stops the script.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
LIST=${1:-$ROOT/meshes.txt}
MESHES=${2:-$HOME/meshes}

[ -f "$LIST" ] || { echo "no mesh file: $LIST"; exit 1; }
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
mkdir -p "$MESHES"

meshValue() { foamDictionary -entry "$2" -value "$1/system/meshDict" 2>/dev/null; }

while read -r NAME REST; do
    case "$NAME" in ''|'#'*) continue ;; esac
    MESH=$MESHES/$NAME

    H=""; FIN=""; NL=""; TR=""
    set -- $REST
    while [ $# -gt 0 ]; do
        case "$1" in
            --maxCellSize)    H=$2;   shift 2 ;;
            --finLevel)       FIN=$2; shift 2 ;;
            --nLayers)        NL=$2;  shift 2 ;;
            --thicknessRatio) TR=$2;  shift 2 ;;
            *)                shift ;;
        esac
    done

    if [ -d "$MESH/constant/polyMesh" ]; then
        if { [ -n "$H" ] && [ "$(meshValue "$MESH" maxCellSize)" != "$H" ]; } ||
           { [ -n "$FIN" ] && [ "$(meshValue "$MESH" localRefinement/fins/additionalRefinementLevels)" != "$FIN" ]; } ||
           { [ -n "$NL" ] && [ "$(meshValue "$MESH" boundaryLayers/patchBoundaryLayers/walls/nLayers)" != "$NL" ]; } ||
           { [ -n "$TR" ] && [ "$(meshValue "$MESH" boundaryLayers/patchBoundaryLayers/walls/thicknessRatio)" != "$TR" ]; }; then
            echo "!! $MESH was built with other options than the line asks, stopping"
            exit 1
        fi
        echo "== $NAME  exists, keeping it"
        continue
    elif [ -e "$MESH" ]; then
        echo "!! $MESH exists but holds no mesh: delete it, stopping"
        exit 1
    fi

    echo "== $NAME  $REST"
    if ! "$ROOT/mesh/Allmesh" "$MESH" $REST < /dev/null; then
        echo "!! $MESH failed to build, stopping"
        exit 1
    fi
done < "$LIST"

echo
echo "done.  meshes under $MESHES"

#!/bin/bash
# ./runCase.sh [--np <N>] [--socket <S>] <sweep-file> [runs-dir] [meshes-dir]
#                                    default: ~/runs, ~/meshes
#
# Makes and runs every case in the sweep file, one after the other.  Each line
# names its mesh with --mesh; a bare name lives in meshes-dir.  The meshes must
# exist: build them first with runMesh.sh.  A case whose directory already
# exists is skipped, so re-running resumes.  --init names the run whose
# solution starts this one; a bare name lives in runs-dir.
#
# --np and --socket apply to every case; an option on a line of the sweep file
# overrides them.  Default --np: the physical cores of one socket.  On a
# computer with several sockets, run one sweep file on each socket.
#
# The script reads the sweep file while it runs: do not edit it in place until
# it ends.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd)
NP=""; SOCKET=""
while [ $# -gt 0 ]; do
    case "$1" in
        --np)     NP=$2; shift 2 ;;
        --socket) SOCKET=$2; shift 2 ;;
        *)        break ;;
    esac
done
SWEEP=${1:?usage: runCase.sh [--np <N>] [--socket <S>] <sweep-file> [runs-dir] [meshes-dir]}
RUNS=${2:-$HOME/runs}
MESHES=${3:-$HOME/meshes}

[ -f "$SWEEP" ] || { echo "no sweep file: $SWEEP"; exit 1; }
[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
if [ -z "$NP" ]; then
    NP=$(lscpu -p=SOCKET,CORE | grep -v '^#' | awk -F, -v s="${SOCKET:-0}" '$1 == s {print $2}' | sort -un | wc -l)
fi
DEFAULTS=(--np "$NP")
[ -n "$SOCKET" ] && DEFAULTS+=(--socket "$SOCKET")
echo "processes for each case: $NP${SOCKET:+, socket $SOCKET}"
mkdir -p "$RUNS"

while read -r NAME REST; do
    case "$NAME" in ''|'#'*) continue ;; esac
    DIR="$RUNS/$NAME"

    if [ -d "$DIR" ]; then
        echo "== $NAME  already exists, skipping"
        continue
    fi

    MESH=""; OPTS=("${DEFAULTS[@]}")
    set -- $REST
    while [ $# -gt 0 ]; do
        case "$1" in
            --mesh) MESH=$2; shift 2 ;;
            --init) I=${2/#\~/$HOME}
                    case "$I" in */*) ;; *) I=$RUNS/$I ;; esac
                    OPTS+=(--init "$I"); shift 2 ;;
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

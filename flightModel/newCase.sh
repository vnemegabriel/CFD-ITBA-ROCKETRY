#!/bin/bash
# ./newCase.sh <run-dir> [options]
#
#     ./newCase.sh ~/runs/m03 --mesh ~/meshes/base --Minf 0.3
#     ./newCase.sh ~/runs/m12 --mesh ~/meshes/base --Minf 1.2 --np 16
#     ./newCase.sh ~/runs/a05 --mesh ~/meshes/base --Minf 0.8 --alpha 5 --refine tip
#
# Copies case-central plus a mesh built by mesh/Allmesh.  See docs/WORKFLOW.md.
# Keep run directories out of OneDrive and out of /mnt/c.
set -e
ROOT=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
[ $# -ge 1 ] || usage
RUN=$1; shift
OPTIONS="$*"

[ -n "$WM_PROJECT_DIR" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
declare -A SET
NP=""; REFINE=""; MESH=""; TEMPLATE=case-central
while [ $# -gt 0 ]; do
    case "$1" in
        --alpha|--beta|--Ti|--nuRatio|--Minf|--pInf|--Tinf)
            SET[${1#--}]=$2; shift 2 ;;
        --mesh)   MESH=$2;   shift 2 ;;
        --np)     NP=$2;     shift 2 ;;
        --refine) REFINE=$2; shift 2 ;;
        *) echo "unknown option $1"; usage ;;
    esac
done

[ -n "$MESH" ] || { echo "--mesh <mesh-dir> is required"; usage; }
[ -d "$MESH/constant/polyMesh" ] || { echo "no mesh in $MESH: run mesh/Allmesh $MESH first"; exit 1; }

case "$REFINE" in ''|tip|fins) ;; *) echo "--refine takes tip or fins"; usage ;; esac

[ -e "$RUN" ] && { echo "!! $RUN exists -- pick a new directory"; exit 1; }
mkdir -p "$RUN"

rsync -a --exclude 'constant/polyMesh' --exclude 'log.*' --exclude '/0' \
      --exclude 'processor*' --exclude 'postProcessing' --exclude '[1-9]*' \
      "$ROOT/$TEMPLATE/" "$RUN/"
MESH=$(cd "$MESH" && pwd)
if [ -n "$REFINE" ]; then
    cp -r "$MESH/constant/polyMesh" "$RUN/constant/"
else
    ln -s "$MESH/constant/polyMesh" "$RUN/constant/polyMesh"
fi
cp "$MESH/meshInfo" "$RUN/constant/meshInfo"
cp "$ROOT/common/Allrun" "$ROOT/common/Allrefine" "$RUN/"
chmod +x "$RUN"/Allrun "$RUN"/Allrefine

cd "$RUN"
touch case.foam
cat > run.info <<EOF
case        "$(basename "$PWD")";
created     "$(date -Iseconds)";
commit      "$(git -C "$ROOT" describe --always --dirty 2>/dev/null || echo unknown)";
mesh        "$MESH";
nCells      $(foamDictionary -entry nCells -value constant/meshInfo);
options     "$OPTIONS";
EOF

TEMPLATES=$(find . -name '*.j2')
if [ -n "$TEMPLATES" ]; then
    command -v jinja2 > /dev/null || { echo "!! needs jinja2:  pipx install jinja2-cli"; exit 1; }
    for f in $TEMPLATES; do
        jinja2 --strict "$f" config.json --format=json -o "${f%.j2}"
        rm "$f"
    done
    echo "rendered $(echo "$TEMPLATES" | wc -l) templates from config.json"
fi

echo "template: $TEMPLATE  ($(foamDictionary -entry application -value system/controlDict))"
for k in "${!SET[@]}"; do
    foamDictionary -entry "$k" -value system/flowConditions > /dev/null 2>&1 || {
        echo "!! $TEMPLATE has no flowConditions entry '$k'."; exit 1; }
    foamDictionary -entry "$k" -set "${SET[$k]}" system/flowConditions > /dev/null
    echo "flowConditions: $k = ${SET[$k]}"
done
[ -n "$NP" ] && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict > /dev/null

[ -n "$REFINE" ] && ./Allrefine "$REFINE"

echo
echo "ready: cd $RUN && ./Allrun ${NP:+$NP}"

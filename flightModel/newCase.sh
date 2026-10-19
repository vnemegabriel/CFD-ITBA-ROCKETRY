#!/bin/bash
# ./newCase.sh <run-dir> [options]
#
#     ./newCase.sh ~/runs/m03 --mesh ~/meshes/base --Minf 0.3
#     ./newCase.sh ~/runs/m12 --mesh ~/meshes/base --Minf 1.2 --np 16
#                  (without --np: every physical core, or those of --socket)
#     ./newCase.sh ~/runs/a05 --mesh ~/meshes/base --Minf 0.8 --alpha 5 --refine tip
#     ./newCase.sh ~/runs/m06 --mesh ~/meshes/base --Minf 0.6 --np 8 --socket 1
#     ./newCase.sh ~/runs/m12 --mesh ~/meshes/base --Minf 1.2 --init ~/runs/m18
#     ./newCase.sh ~/runs/m03 --mesh ~/meshes/base --Minf 0.3 --solver pimple
#     ./newCase.sh ~/runs/f12 --mesh ~/meshes/h0 --Minf 1.2 --solver forced --init ~/runs/m12
#                  --flutterModes ~/flutter/pulse1   (fin_flutter: python -m finflutter modes)
#
# Copies case-<solver> (default central) plus a mesh built by mesh/Allmesh.
# A template other than central holds only what differs from its parent
# (file parent, default central); <parent>.exclude lists the parent's files
# it drops.
# See docs/guide, chapter Cases.
# Keep run directories out of /mnt/c.
set -e
ROOT=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
[ $# -ge 1 ] || usage
RUN=$1; shift
OPTIONS="$*"

[ -n "$WM_PROJECT_DIR" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
declare -A SET
NP=""; SOCKET=""; REFINE=""; MESH=""; INIT=""; SOLVER=central; MODES=""
while [ $# -gt 0 ]; do
    case "$1" in
        --alpha|--beta|--Ti|--nuRatio|--Minf|--pInf|--Tinf)
            SET[${1#--}]=$2; shift 2 ;;
        --mesh)   MESH=$2;   shift 2 ;;
        --np)     NP=$2;     shift 2 ;;
        --socket) SOCKET=$2; shift 2 ;;
        --refine) REFINE=$2; shift 2 ;;
        --init)   INIT=$2;   shift 2 ;;
        --solver) SOLVER=$2; shift 2 ;;
        --flutterModes) MODES=${2/#\~/$HOME}; shift 2 ;;
        *) echo "unknown option $1"; usage ;;
    esac
done

[ -n "$MESH" ] || { echo "--mesh <mesh-dir> is required"; usage; }
[ -d "$MESH/constant/polyMesh" ] || { echo "no mesh in $MESH: run mesh/Allmesh $MESH first"; exit 1; }

case "$REFINE" in ''|tip|fins) ;; *) echo "--refine takes tip or fins"; usage ;; esac
[ -d "$ROOT/case-$SOLVER" ] || { echo "--solver: no template $ROOT/case-$SOLVER"; usage; }
TEMPLATE=case-$SOLVER
if [ -n "$INIT" ]; then
    [ -d "$INIT/system" ] || { echo "--init $INIT is not a run directory"; exit 1; }
    INIT=$(cd "$INIT" && pwd)
fi

[ -e "$RUN" ] && { echo "!! $RUN exists -- pick a new directory"; exit 1; }
mkdir -p "$RUN"

COPY=(rsync -a --exclude 'constant/polyMesh' --exclude 'log.*' --exclude '/0'
      --exclude 'processor*' --exclude 'postProcessing' --exclude '[1-9]*')
chain() {
    if [ "$1" = central ]; then echo central; return; fi
    chain "$(cat "$ROOT/case-$1/parent" 2>/dev/null || echo central)"
    echo "$1"
}
PARENT=""
for T in $(chain "$SOLVER"); do
    if [ -n "$PARENT" ]; then
        while read -r X; do
            [ -n "$X" ] && rm -rf "${RUN:?}/${X#/}"
        done < "$ROOT/case-$T/$PARENT.exclude"
    fi
    "${COPY[@]}" --exclude '/*.exclude' --exclude '/parent' "$ROOT/case-$T/" "$RUN/"
    PARENT=$T
done
MESH=$(cd "$MESH" && pwd)
if [ -n "$REFINE" ]; then
    cp -r "$MESH/constant/polyMesh" "$RUN/constant/"
else
    ln -s "$MESH/constant/polyMesh" "$RUN/constant/polyMesh"
fi
cp "$MESH/meshInfo" "$RUN/constant/meshInfo"
if [ -n "$MODES" ]; then
    [ -f "$MODES" ] || { echo "!! no such flutterModes file: $MODES"; exit 1; }
    cp "$MODES" "$RUN/constant/flutterModes"
    (cd "$ROOT/../fin_flutter" && pwd) > "$RUN/finflutterPath"
elif [ -f "$RUN/system/flutterModes.H" ]; then
    echo "!! $TEMPLATE needs constant/flutterModes: --flutterModes <file> (python -m finflutter modes)"; exit 1
fi
cp "$ROOT/common/Allrun" "$ROOT/common/Allrefine" "$ROOT/common/machineLog.py" "$RUN/"
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
[ -n "$INIT" ] && echo "$INIT" > initFrom && echo "init        \"$INIT\";" >> run.info

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
# default: every physical core of the computer, or of the socket of --socket
if [ -z "$NP" ]; then
    NP=$(lscpu -p=SOCKET,CORE | grep -v '^#' | awk -F, -v s="$SOCKET" 's == "" || s == "none" || $1 == s {print $1","$2}' | sort -u | wc -l)
fi
foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict > /dev/null
echo "processes: $NP"
[ -n "$SOCKET" ] && echo "$SOCKET" > socket

[ -n "$REFINE" ] && ./Allrefine "$REFINE"

meshHash() { (cd "$1" && cat points faces owner neighbour boundary | sha256sum | cut -c1-16); }
if [ -n "$REFINE" ]; then
    HASH=$(meshHash constant/polyMesh)
elif [ "$MESH/meshHash" -nt "$MESH/constant/polyMesh/points" ]; then
    HASH=$(cat "$MESH/meshHash")
else
    HASH=$(meshHash "$MESH/constant/polyMesh")
    { echo "$HASH" > "$MESH/meshHash.$$" && mv "$MESH/meshHash.$$" "$MESH/meshHash"; } 2> /dev/null || true
fi
TH=constant/thermophysicalProperties
MW=$(foamDictionary -entry mixture/specie/molWeight -value $TH)
CP=$(foamDictionary -entry mixture/thermodynamics/Cp -value $TH)
fd() { foamDictionary -entry "$1" -value "$2" 2>/dev/null || echo none; }
{
    echo "meshHash    \"$HASH\";"
    echo "flow"
    echo "{"
    foamDictionary -expand system/flowDerived | sed -n '/^}/,$p' | sed '1d; /^$/d; /^\/\//d; s/^/    /'
    LC_ALL=C awk -v mw="$MW" -v cp="$CP" -v r="$(fd rhoInf system/flowDerived)" -v u="$(fd Uinf system/flowDerived)" \
        -v mu="$(fd muInf system/flowDerived)" -v d="$(fd lRef constant/meshInfo)" \
        -v l="$(fd lBody constant/meshInfo)" 'BEGIN{
        R = 8314.46/mw; cv = cp - R
        printf "    RgasThermo      %.6g;\n", R
        printf "    PrEucken        %.4f;\n", cp/(cv*(1.32 + 1.77*R/cv))
        printf "    RePerMetre      %.6g;\n", r*u/mu
        printf "    ReD             %.6g;\n", r*u*d/mu
        printf "    ReL             %.6g;\n", r*u*l/mu }'
    echo "    thermo          \"$(fd thermoType/type $TH) $(fd thermoType/transport $TH) $(fd thermoType/energy $TH)\";"
    echo "    sutherland      ($(fd mixture/transport/As $TH) $(fd mixture/transport/Ts $TH));"
    echo "    Cp              $CP;"
    echo "    turbulence      $(fd RAS/RASModel constant/turbulenceProperties);"
    echo "    wallT           $(fd boundaryField/walls/type 0.orig/T);"
    echo "    wallNut         $(fd boundaryField/walls/type 0.orig/nut);"
    echo "    initial         \"${INIT:-uniform free stream}\";"
    echo "}"
} >> run.info

echo
echo "ready: cd $RUN && ./Allrun ${NP:+$NP}"

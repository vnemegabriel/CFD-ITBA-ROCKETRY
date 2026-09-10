#!/bin/bash
# ./newCase.sh <run-dir> [options]
#
#     ./newCase.sh ~/runs/m02 --regime sub   --Uinf 68
#     ./newCase.sh ~/runs/m12 --regime super --Minf 1.2 --np 16
#     ./newCase.sh ~/runs/a05 --regime sub   --Uinf 68 --alpha 5 --refine tip
#
# Copies a template plus the mesh built by mesh/Allmesh.  See docs/WORKFLOW.md.
# Keep run directories out of OneDrive and out of /mnt/c.
set -e
ROOT=$(cd "$(dirname "$0")" && pwd)

usage() { sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
[ $# -ge 1 ] || usage
RUN=$1; shift

[ -n "$WM_PROJECT_DIR" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
[ -d "$ROOT/mesh/constant/polyMesh" ] || { echo "no mesh: run mesh/Allmesh first"; exit 1; }

declare -A SET
NP=""; REFINE=""; REGIME="sub"
while [ $# -gt 0 ]; do
    case "$1" in
        --alpha|--beta|--Uinf|--nu|--rhoInf|--Ti|--nuRatio|--Minf|--pInf|--Tinf)
            SET[${1#--}]=$2; shift 2 ;;
        --np)     NP=$2;     shift 2 ;;
        --refine) REFINE=$2; shift 2 ;;
        --regime) REGIME=$2; shift 2 ;;
        *) echo "unknown option $1"; usage ;;
    esac
done

case "$REFINE" in ''|tip|fins) ;; *) echo "--refine takes tip or fins"; usage ;; esac
case "$REGIME" in
    sub)   TEMPLATE=case-subsonic ;;
    trans) TEMPLATE=case-transonic ;;
    super) TEMPLATE=case-supersonic ;;
    *) echo "--regime takes sub, trans or super"; usage ;;
esac

[ -e "$RUN" ] && { echo "!! $RUN exists -- pick a new directory"; exit 1; }
mkdir -p "$RUN"

rsync -a --exclude 'constant/polyMesh' --exclude 'log.*' --exclude '/0' \
      --exclude 'processor*' --exclude 'postProcessing' --exclude '[1-9]*' \
      "$ROOT/$TEMPLATE/" "$RUN/"
cp -r "$ROOT/mesh/constant/polyMesh" "$RUN/constant/"
cp "$ROOT/mesh/meshInfo" "$RUN/constant/meshInfo"
cp "$ROOT/common/Allrun" "$ROOT/common/Allrefine" "$RUN/"
chmod +x "$RUN"/Allrun "$RUN"/Allrefine

cd "$RUN"
touch case.foam

echo "template: $TEMPLATE  ($(foamDictionary -entry application -value system/controlDict))"
for k in "${!SET[@]}"; do
    foamDictionary -entry "$k" -value system/flowConditions > /dev/null 2>&1 || {
        echo "!! $TEMPLATE has no flowConditions entry '$k'."
        echo "   sub takes Uinf/nu/rhoInf; trans and super take Minf/pInf/Tinf."
        exit 1; }
    foamDictionary -entry "$k" -set "${SET[$k]}" system/flowConditions > /dev/null
    echo "flowConditions: $k = ${SET[$k]}"
done
[ -n "$NP" ] && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict > /dev/null

[ -n "$REFINE" ] && ./Allrefine "$REFINE"

echo
echo "ready: cd $RUN && ./Allrun ${NP:+$NP}"

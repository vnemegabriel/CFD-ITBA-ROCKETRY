#!/bin/bash
# newCase.sh -- start a run from the case template, in one line.
#
#     ./newCase.sh ~/runs/a05  mesh/output/aconcagua_half_s2.msh  --alpha 5
#     ./newCase.sh ~/runs/ref  mesh/output/aconcagua_quarter_s1.msh
#     ./newCase.sh ~/runs/roll mesh/output/aconcagua_full_s2.msh  --alpha 4 --beta 4 --np 16
#     ./newCase.sh ~/runs/tip  mesh/output/aconcagua_half_s2.msh  --alpha 5 --refine tip
#     ./newCase.sh ~/runs/m08  mesh/output/aconcagua_half_s2.msh  --regime trans --Minf 0.8
#
# Copies ONE OF THE TEMPLATES (never run in place) to the run directory, sets
# the flow conditions you pass, converts the mesh, optionally refines around
# the fins, and tells you what to run.  Keep run directories OUT of OneDrive
# and, under WSL, out of /mnt/c: both make OpenFOAM's I/O several times slower.
#
# Options:
#   --regime sub|trans|super   which template.  Default sub.
#                                sub    case-subsonic    simpleFoam       M < 0.3
#                                trans  case-transonic   rhoSimpleFoam    0.3-1.2
#                                super  case-supersonic  rhoCentralFoam   M > 1.2
#   --np N                     decomposition
#   --refine tip|fins          local refinement at the fins, see Allrefine
#
#   flow conditions, written into system/flowConditions:
#     every template   --alpha --beta --Ti --nuRatio
#     sub only         --Uinf --nu --rhoInf
#     trans/super only --Minf --pInf --Tinf
#
# Needs OpenFOAM sourced:  source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -e
ROOT=$(cd "$(dirname "$0")" && pwd)

usage() {
    sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
}
[ $# -ge 2 ] || usage
RUN=$1; MSH=$2; shift 2
[ -f "$MSH" ] || { echo "no such mesh: $MSH"; exit 1; }
[ -n "$WM_PROJECT_DIR" ] || { echo "source the OpenFOAM bashrc first"; exit 1; }
MSH=$(cd "$(dirname "$MSH")" && pwd)/$(basename "$MSH")

declare -A SET
NP=""
REFINE=""
REGIME="sub"
while [ $# -gt 0 ]; do
    case "$1" in
        --alpha|--beta|--Uinf|--nu|--rhoInf|--Ti|--nuRatio|--Minf|--pInf|--Tinf)
            SET[${1#--}]=$2; shift 2 ;;
        --np) NP=$2; shift 2 ;;
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
[ -d "$ROOT/$TEMPLATE" ] || { echo "no template $ROOT/$TEMPLATE"; exit 1; }

if [ -e "$RUN" ]; then
    echo "!! $RUN exists -- pick a new directory (this script never overwrites a run)"
    exit 1
fi
mkdir -p "$RUN"
# the template only: no mesh, fields, logs or results
rsync -a --exclude 'constant/polyMesh' --exclude 'log.*' --exclude '/0' \
      --exclude 'processor*' --exclude 'postProcessing' --exclude '[1-9]*' \
      "$ROOT/$TEMPLATE/" "$RUN/"
cd "$RUN"

echo "template: $TEMPLATE  ($(foamDictionary -entry application -value system/controlDict))"
for k in "${!SET[@]}"; do
    # refuse silently-ignored settings: e.g. --Uinf on a compressible template,
    # where the velocity is DERIVED from Minf and the ambient state
    foamDictionary -entry "$k" -value system/flowConditions > /dev/null 2>&1 || {
        echo "!! $TEMPLATE/system/flowConditions has no entry '$k'."
        echo "   sub takes Uinf/nu/rhoInf; trans and super take Minf/pInf/Tinf."
        exit 1; }
    foamDictionary -entry "$k" -set "${SET[$k]}" system/flowConditions > /dev/null
    echo "flowConditions: $k = ${SET[$k]}"
done
[ -n "$NP" ] && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict > /dev/null

./Allmesh "$MSH"
if [ -n "$REFINE" ]; then ./Allrefine "$REFINE"; fi

echo
echo "run directory: $RUN"
echo "   cd $RUN && ./Allrun ${NP}"

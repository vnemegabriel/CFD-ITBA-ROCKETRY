#!/bin/bash
# newCase.sh -- start a run from the case template, in one line.
#
#     ./newCase.sh ~/runs/a05  mesh/output/aconcagua_half_s2.msh  --alpha 5
#     ./newCase.sh ~/runs/ref  mesh/output/aconcagua_quarter_s1.msh
#     ./newCase.sh ~/runs/roll mesh/output/aconcagua_full_s2.msh  --alpha 4 --beta 4 --np 16
#     ./newCase.sh ~/runs/tip  mesh/output/aconcagua_half_s2.msh  --alpha 5 --refine tip
#
# Copies case/ (the template, never run in place) to the run directory, sets
# the flow conditions you pass, converts the mesh, optionally refines around
# the fins, and tells you what to run.  Keep run directories OUT of OneDrive
# and, under WSL, out of /mnt/c: both make OpenFOAM's I/O several times slower.
#
# Options:  --alpha --beta --Uinf --nu --rhoInf --Ti --nuRatio   flow conditions
#           --np N                                               decomposition
#           --refine tip|fins                                    see case/Allrefine
#
# Needs OpenFOAM sourced:  source /usr/lib/openfoam/openfoam2412/etc/bashrc
set -e
ROOT=$(cd "$(dirname "$0")" && pwd)

usage() {
    sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
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
while [ $# -gt 0 ]; do
    case "$1" in
        --alpha|--beta|--Uinf|--nu|--rhoInf|--Ti|--nuRatio) SET[${1#--}]=$2; shift 2 ;;
        --np) NP=$2; shift 2 ;;
        --refine) REFINE=$2; shift 2 ;;
        *) echo "unknown option $1"; usage ;;
    esac
done
case "$REFINE" in ''|tip|fins) ;; *) echo "--refine takes tip or fins"; usage ;; esac

if [ -e "$RUN" ]; then
    echo "!! $RUN exists -- pick a new directory (this script never overwrites a run)"
    exit 1
fi
mkdir -p "$RUN"
# the template only: no mesh, fields, logs or results
rsync -a --exclude 'constant/polyMesh' --exclude 'log.*' --exclude '/0' \
      --exclude 'processor*' --exclude 'postProcessing' --exclude '[1-9]*' \
      "$ROOT/case/" "$RUN/"
cd "$RUN"

for k in "${!SET[@]}"; do
    foamDictionary -entry "$k" -set "${SET[$k]}" system/flowConditions > /dev/null
    echo "flowConditions: $k = ${SET[$k]}"
done
[ -n "$NP" ] && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict > /dev/null

./Allmesh "$MSH"
if [ -n "$REFINE" ]; then ./Allrefine "$REFINE"; fi

echo
echo "run directory: $RUN"
echo "   cd $RUN && ./Allrun ${NP}"

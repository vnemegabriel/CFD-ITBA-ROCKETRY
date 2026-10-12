#!/bin/sh
# ./post/sweepTable.sh [runs-dir] > sweep.csv          default: ~/runs
#
# One CSV line per run in runs-dir: flow conditions, the last run.info
# record, the coefficients of the last line of forceCoeffs1, and y+ on the
# walls at the last time.  A run that has not started has empty fields.
RUNS=${1:-$HOME/runs}

value() { awk -v k="$1" '$1 == k { v = $2 } END { gsub(/[";]/, "", v); print v }' "$2" 2>/dev/null; }
last() { d=$(ls "$1" 2>/dev/null | sort -g | tail -n 1); [ -n "$d" ] && echo "$1/$d/"; }

echo "case,mesh,nCells,Minf,pInf,Tinf,alpha,status,steps,tEnd,wallSeconds,Cd,Cl,CmPitch,yPlusMin,yPlusMax"
for R in "$RUNS"/*/; do
    R=${R%/}
    [ -f "$R/run.info" ] || continue
    FC=$R/system/flowConditions
    INFO=$R/run.info
    MESH=$(basename "$(value mesh "$INFO")")
    TEND=$(awk '$1 == "simTime" { gsub(/[();]/, ""); t = $3 } END { print t }' "$INFO")

    C=$(last "$R/postProcessing/forceCoeffs1")coefficient.dat
    COEF=$(awk '/^# Time/ { for (i = 2; i <= NF; i++) c[$i] = i - 1 } !/^#/ { l = $0 }
        END { if (l != "") { split(l, f); printf "%s,%s,%s", f[c["Cd"]], f[c["Cl"]], f[c["CmPitch"]] } else printf ",," }' "$C" 2>/dev/null || echo ",,")

    Y=$(last "$R/postProcessing/yPlus")yPlus.dat
    YP=$(awk '!/^#/ { t = $1; if (!(t in mn) || $3 < mn[t]) mn[t] = $3; if (!(t in mx) || $4 > mx[t]) mx[t] = $4 }
        END { if (t != "") printf "%s,%s", mn[t], mx[t]; else printf "," }' "$Y" 2>/dev/null || echo ",")

    echo "$(basename "$R"),$MESH,$(value nCells "$INFO"),$(value Minf "$FC"),$(value pInf "$FC"),$(value Tinf "$FC"),$(value alpha "$FC"),$(value status "$INFO"),$(value steps "$INFO"),$TEND,$(value wallSeconds "$INFO"),$COEF,$YP"
done

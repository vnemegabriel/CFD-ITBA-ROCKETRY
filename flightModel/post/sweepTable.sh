#!/bin/sh
# ./post/sweepTable.sh [runs-dir] > sweep.csv          default: ~/runs
#
# One CSV line per run in runs-dir: solver, flow conditions, the last run.info
# record, the coefficients of the last line of forceCoeffs1 and their mean over
# the stop-criterion window, the last initial residuals, and y+ on the walls at
# the last time.  A run that has not started has empty fields.
LC_ALL=C; export LC_ALL
RUNS=${1:-$HOME/runs}

value() { awk -v k="$1" '$1 == k { v = $2 } END { gsub(/[";]/, "", v); print v }' "$2" 2>/dev/null; }
last() { d=$(ls "$1" 2>/dev/null | sort -g | tail -n 1); [ -n "$d" ] && echo "$1/$d/"; }

echo "case,mesh,nCells,solver,Minf,pInf,Tinf,alpha,status,steps,tEnd,wallSeconds,Cd,Cl,CmPitch,window,CdMean,ClMean,CmPitchMean,resUx,resUy,resUz,resE,resP,resK,resOmega,yPlusMin,yPlusMax,yPlusAvgWalls,yPlusAvgFins"
for R in "$RUNS"/*/; do
    R=${R%/}
    [ -f "$R/run.info" ] || continue
    FC=$R/system/flowConditions
    INFO=$R/run.info
    MESH=$(basename "$(value mesh "$INFO")")
    TEND=$(awk '$1 == "simTime" { gsub(/[();]/, ""); t = $3 } END { print t }' "$INFO")
    SOLVER=$(cat "$R/system/controlDict" "$R/system/controlDictBase" 2>/dev/null | value application /dev/stdin)
    WIN=$(tr -d ' \n' < "$R/config.json" 2>/dev/null | sed -n 's/.*"coeffs_window":\[\([0-9.eE+-]*\).*/\1/p')

    C=$(last "$R/postProcessing/forceCoeffs1")coefficient.dat
    COEF=$(awk -v w="$WIN" '/^# Time/ { for (i = 2; i <= NF; i++) c[$i] = i - 1; next }
        !/^#/ && NF { n++; t[n] = $1; d[n] = $c["Cd"]; l[n] = $c["Cl"]; m[n] = $c["CmPitch"] }
        END {
            if (!n) { printf ",,,,,,"; exit }
            k = 0; sd = sl = sm = 0
            for (i = n; i >= 1 && (w == "" || t[n] - t[i] <= w); i--) { k++; sd += d[i]; sl += l[i]; sm += m[i] }
            printf "%s,%s,%s,%s,%.6g,%.6g,%.6g", d[n], l[n], m[n], w, sd/k, sl/k, sm/k
        }' "$C" 2>/dev/null || echo ",,,,,,")

    S=$(last "$R/postProcessing/residuals")solverInfo.dat
    RES=$(awk '/^# Time/ { for (i = 2; i <= NF; i++) c[$i] = i - 1; next } !/^#/ && NF { l = $0 }
        END { split(l, f); n = split("Ux Uy Uz e p k omega", v, " ")
              for (i = 1; i <= n; i++) printf "%s%s", (i > 1 ? "," : ""), ((v[i] "_initial") in c ? f[c[v[i] "_initial"]] : "") }' "$S" 2>/dev/null || echo ",,,,,,")

    Y=$(last "$R/postProcessing/yPlus")yPlus.dat
    YP=$(awk '!/^#/ && NF { if ($1 != t) { t = $1; mn = mx = aw = af = "" }
                     if (mn == "" || $3 < mn) mn = $3; if (mx == "" || $4 > mx) mx = $4
                     if ($2 == "walls") aw = $5; if ($2 == "fins") af = $5 }
        END { printf "%s,%s,%s,%s", mn, mx, aw, af }' "$Y" 2>/dev/null || echo ",,,")

    echo "$(basename "$R"),$MESH,$(value nCells "$INFO"),$SOLVER,$(value Minf "$FC"),$(value pInf "$FC"),$(value Tinf "$FC"),$(value alpha "$FC"),$(value status "$INFO"),$(value steps "$INFO"),$TEND,$(value wallSeconds "$INFO"),$COEF,$RES,$YP"
done

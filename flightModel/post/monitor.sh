#!/bin/sh
# ./post/monitor.sh <run-dir> [seconds]
#
# Cd and Cl against the iteration, and the initial residuals, from postProcessing/ of a
# running case; redrawn every [seconds] (default 30).  Joins the directories
# of a continued run.  Without a display it writes <run-dir>/monitor.png.
# Needs gnuplot.
LC_ALL=C; export LC_ALL
RUN=${1:?usage: monitor.sh <run-dir> [seconds]}
EVERY=${2:-30}
PP=$RUN/postProcessing
command -v gnuplot > /dev/null || { echo "needs gnuplot:  sudo apt install gnuplot"; exit 1; }
[ -d "$PP/forceCoeffs1" ] || { echo "no $PP/forceCoeffs1 yet"; exit 1; }

# data column of each name in the '# Time ...' header line
cols() {
    awk -v want="$2" '/^# Time/ { for (i = 2; i <= NF; i++) if ($i ~ want) printf "%d ", i - 1; exit }' "$1"
}
first() { ls "$PP/$1"/*/"$2" 2>/dev/null | head -n 1; }

COEF=$(first forceCoeffs1 coefficient.dat)
RES=$(first residuals solverInfo.dat)
CD=$(cols "$COEF" '^Cd$'); CL=$(cols "$COEF" '^Cl$')
RCOLS=$(cols "$RES" '_initial$')
RNAMES=$(awk '/^# Time/ { for (i = 2; i <= NF; i++) if ($i ~ /_initial$/) { sub(/_initial/, "", $i); printf "%s ", $i }; exit }' "$RES")

if [ -n "$DISPLAY" ]; then
    TERM_SET="set term qt noraise"
else
    TERM_SET="set term pngcairo size 1000,800; set output '$RUN/monitor.png'"
    echo "no display: writing $RUN/monitor.png every $EVERY s"
fi

gnuplot -persist -e "
coef = '< cat $PP/forceCoeffs1/*/coefficient.dat | grep -v \"^#\" | sort -g';
res  = '< cat $PP/residuals/*/solverInfo.dat | grep -v \"^#\" | sort -g';
rcols = '$RCOLS'; rnames = '$RNAMES';
$TERM_SET;
while (1) {
    set multiplot layout 2,1 title '$(basename "$RUN")';
    unset logscale y; set xlabel 'iteration'; set ylabel 'coefficient'; set grid;
    plot coef using 1:$CD with lines title 'Cd', coef using 1:$CL with lines title 'Cl';
    set logscale y; set ylabel 'initial residual';
    plot for [i = 1:words(rcols)] res using 1:(column(int(word(rcols, i)))) with lines title word(rnames, i);
    unset multiplot;
    pause $EVERY;
    if (strlen(GPVAL_TERM) > 0 && GPVAL_TERM eq 'pngcairo') { set output '$RUN/monitor.png' }
}"

mkdir outplot/zerocut

root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
ZeroPlot()
.q
EOF

mkdir outplot/invcut

root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
ZeroProc()
.q
EOF

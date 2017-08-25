mkdir outplot/trkrcorr
mkdir outplot/ecalcorr

root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
TrkrRadCorr()
EcalRadCorr()
.q
EOF

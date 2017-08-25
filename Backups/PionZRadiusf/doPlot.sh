#!/bin/bash

mkdir -p outplot

root -b -l <<EOF
.L style.cxx+
.L RadiusPlot.C+
ProcCmp()
ZeroPlot()
ERCorr()
HeadPlot()
.q
EOF

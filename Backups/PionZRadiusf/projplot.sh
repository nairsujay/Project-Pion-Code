#!/bin/bash

mkdir -p outplot/tailplot
mkdir -p outplot/headplot
mkdir -p outplot/allheadplot
mkdir -p outplot/alltailplot

root -b -l <<EOF
.L style.cxx+
.L RadiusPlot.C+
ProjPlot("Proc2DPlots.root","tailplot");
ProjPlot("Proc2DHeadPlots.root","headplot");
AllProjPlot("Proc2DPlots.root","alltailplot");
AllProjPlot("Proc2DHeadPlots.root","allheadplot");
.q
EOF


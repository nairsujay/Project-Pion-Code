#!/bin/bash

mkdir outplot/trkrvsecal

root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
TrkrvsEcal()
.q
EOF

#!/bin/bash

root -b -l <<EOF
.L style.cxx+
.L RadiusPlot.C+
CutPlot()
.q
EOF

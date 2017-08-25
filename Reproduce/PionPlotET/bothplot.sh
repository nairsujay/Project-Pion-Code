#!/bin/bash

root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
BothPlot()
.q
EOF

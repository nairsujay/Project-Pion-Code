#!/bin/bash

mkdir outplot/meanplots

root -b -l <<EOF
.L style.cxx+
.L Mean.C+
ProjClus()
.q
EOF


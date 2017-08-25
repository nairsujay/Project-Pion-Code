#!/bin/bash


mkdir -p outplot/proccorr

root -b -l <<EOF
TFile::Open("$1")
.x procCorr.C
.q
EOF


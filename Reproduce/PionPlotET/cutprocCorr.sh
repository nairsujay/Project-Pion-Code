#!/bin/bash


mkdir -p outplot/cutproccorr

root -b -l <<EOF
TFile::Open("$1")
.x cutprocCorr.C
.q
EOF


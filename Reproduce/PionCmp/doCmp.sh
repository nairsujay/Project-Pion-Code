#!/bin/bash


mkdir outplot

root -b -l <<EOF
TFile::Open("ZEnergyDepositbyProcess.root")
.x doCmp2.C
.q
EOF


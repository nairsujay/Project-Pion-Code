#!/bin/bash

mkdir -p outplot/reproduce
mkdir -p outplot/zerocut 

root -b -l <<EOF
.L style.cxx+ 
.L ReproducePlot.C+
ReproducePlot()
.q
EOF


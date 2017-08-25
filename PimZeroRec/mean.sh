#!/bin/bash

root -b -l <<EOF
.L style.cxx+
.L Mean.C+
BasicMean();
RadiusMean();
.q
EOF

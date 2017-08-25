#!/bin/bash

root -b -l <<EOF
.L style.cxx+
.L RadiusMean.C+
RadiusMean()
.q
EOF


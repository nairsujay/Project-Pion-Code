#!/bin/bash

rm -f Pion_Output_Recur.root
cp PionM_Output_MINERvA.root Pion_Output_Recur.root

root -b -l <<EOF
.L ZeroRec.C+
ZeroRec()
.q
EOF

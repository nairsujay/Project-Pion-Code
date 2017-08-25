#!/bin/bash
root -b -l <<EOF
.L dEdxAnaTool.C+
dEdxAnaTool Tes
Tes.Loop("$1")
.q
EOF

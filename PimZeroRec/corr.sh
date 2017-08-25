#!/bin/bash

mkdir outplot/ERcorrplots
mkdir outplot/ZTcorrplots
mkdir outplot/ZEcorrplots
mkdir outplot/TEcorrplots
mkdir outplot/NEcorrplots

root -b -l <<EOF
.L style.cxx+
.L ERplot.C+
CorrPlot("Radiusunbin","Energyunbin","Radius against Energy;Radius (cm);Energy (MeV)",20,0,100,60,0,30,"ERcorrplots",kTRUE)
CorrPlot("Zarray","Tarray","Z against T;Z (mm);T (mm)",100,5000,10000,44,-1100,1100,"ZTcorrplots",kTRUE)
CorrPlot("Zarray","Energyunbin","Z against Energy;Z (mm);Energy (MeV)",100,5000,10000,60,0,30,"ZEcorrplots",kTRUE)
CorrPlot("Tarray","Energyunbin","T against Energy;T (mm);Energy (MeV)",44,-1100,1100,60,0,30,"TEcorrplots",kTRUE)
CorrPlot("neighrad","Energyunbin","Distance to Neighbour against Energy;Distance (cm);Energy(MeV)",25,0,50,60,0,30,"NEcorrplots",kTRUE)
.q
EOF

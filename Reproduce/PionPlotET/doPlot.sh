#!/bin/bash

mkdir outplot/reproduce
mkdir outplot/cutproccorr
mkdir outplot/multicut

#Possible Cuts
singlecut="(trkr30<25)"
senglecut="(trkr30<20)"
peakcut="(trkr30<50)&&(trkr25<40)&&(trkr20<30)&&(trkr15<20)&&(trkr10<10)" 
strictcut="(trkr30<40)&&(trkr25<30)&&(trkr20<20)&&(trkr15<10)&&(trkr10<10)"
flatcut="(trkr30<20)&&(trkr25<20)&&(trkr20<20)&&(trkr15<20)&&(trkr10<20)"
strictflatcut="(trkr30<05)&&(trkr25<05)&&(trkr20<05)&&(trkr15<05)&&(trkr10<05)"
abscut="(trkr30<05)&&(trkr25<05)&&(trkr20<03)&&(trkr15<02)&&(trkr10<01)"
vstrictcut="(trkr30<01)&&(trkr25<.5)&&(trkr20<.4)&&(trkr15<.3)&&(trkr10<.1)"
corrcut="(trkr30<01)&&(trkr25<.5)&&(trkr20<.4)&&(trkr15<.3)&&(trkr10<.1)&&(trkr10==trkr15)"
vvstrictcut="(trkr30<.1)&&(trkr25<.05)&&(trkr20<.04)&&(trkr15<.03)&&(trkr10<.01)"
trialcut="(trkr30<25)&&(trkr10>10)"


#So far after zero suppression 
#(trkr30<25): rangeout->~25% [Suspect this is optimum (confim after correct zero suppression)]
#(trkr30<20): rangeout->~25% [Worse than above]
#(trkr30<25)&&(trkr10<20): rangeout->~25% [Worse than above (increases abs)]
#(trkr30<20)&&(trkr10>10): [Attempt to decrease abs] rangeout->~30% [some success but increases inel ------OPTIMIZE THIS]       ----HEAVILY DEPENDENT ON PHYSICAL ZEROS and SUPPRESSION





if [ $1 -eq 1 ]
then
root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
ReproducePlot()
.q
EOF
fi


if [ $2 -eq 1 ]
then
root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
CutPlot(50,20,30)
.q
EOF
fi

if [ $3 -eq 1 ]
then
root -b -l <<EOF
.L style.cxx+
.L doPlot.C+
MultiCut("$singlecut","$trialcut")
.q
EOF
fi


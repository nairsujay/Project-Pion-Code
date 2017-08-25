#include "style.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "THStack.h"
#include "PlotHelp.h"


void ReproducePlot()
{
    TList* lout = new TList ;

    const TString Proc[8] = {"Absorption","Multi-Production","Charge Exchange","Double Charge Exchange","Inelastic","Rangeout","Decay","all"};
    const TString zProc[8] = {"zAbsorption","zMulti-Production","zCharge Exchange","zDouble Charge Exchange","zInelastic","zRangeout","zDecay","zall"};
    const TString Cut[8] = {"(secondary_proc==0)","(secondary_proc==1)","(secondary_proc==2)",
        "(secondary_proc==3)","(secondary_proc==4)","(secondary_proc==5)",
        "(secondary_proc==6)","(pdg==-211)"};
    // pdg cut does not remove any data, 
    //it's purely for syntax as Plot1DModule requires a cut
    //
    const TString zerocut("&&(trkr + ecal > 0)");
    const TString ozcut("&&(trkr10 + ecal10 > 0)");

    const TString llrcut("&&(llrproton < 0)");

    TTree* tree = PlotHelp::InputFile("PionM_Output_MINERvA.root","tout");

    TCanvas* c1 = new TCanvas("c1","",600,400);

    for(Int_t ii=0;ii<8;ii++) //Nested loops iterate over process and radius
    {
        for(Int_t jj=10;jj<=30;jj+=5)
        {        
            PlotHelp::Plot1DModule(tree,c1,Proc[ii]+Form("_%d",jj),Form("trkr%d/0.8+ecal%d/0.45",jj,jj),Cut[ii]+llrcut,30,0,150,lout,"outplot/reproduce",ii);
            PlotHelp::Plot1DModule(tree,c1,zProc[ii]+Form("_%d",jj),Form("trkr%d/0.8+ecal%d/0.45",jj,jj),Cut[ii]+llrcut+Form("&&((trkr%d/0.8+ecal%d/0.45)>0)",jj,jj),30,0,150,lout,"outplot/reproduce",ii);

        }
    }

    TFile* fout =new TFile("outplot/reproduce/EnergyDepositbyProcess.root","RECREATE");
    lout->Write();
    fout->Close();

    delete lout;
    delete fout;
    delete c1;
}

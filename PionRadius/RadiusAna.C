#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2D.h"
#include "TList.h"

#include "PlotHelp.h"

#include <cmath>


TTree* InputFile(const TString fname,const TString tname)
{
    TFile* fin = new TFile(fname,"READ");
    TTree* tree =(TTree*)fin->Get(tname);
    return tree;
}



void Printer()
{
    TTree* tree = InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");
    Int_t radiusrange;
    Int_t Radiusbin[2000];
    tree->SetBranchAddress("radiusrange",&radiusrange);
    tree->SetBranchAddress("Radiusbin",&Radiusbin);

    for(Int_t ii = 0; ii<tree->GetEntries(); ii++)
    {
        tree->GetEntry(ii);
        if(Radiusbin[0] == 0){printf("\n\nFailed\n\n");exit(1);}; //exit if filtering failed
        for(Int_t jj =0; jj < radiusrange; jj++){printf(" %d ",Radiusbin[jj]);};
        printf("\n");
    }

    delete tree;

}

void ZeroFilter()
{
    const TString trkrcut("&&(ecal==0)&&(ecal10==0)&&(ecal15==0)&&(ecal20==0)&&(ecal25==0)&&(ecal30==0)");
    TTree* intree = InputFile("Pion_Output_MINERvA.root","tout");
    Int_t radiusrange;
    Int_t Radiusbin[2000];
    intree->SetBranchAddress("radiusrange",&radiusrange);
    intree->SetBranchAddress("Radiusbin",&Radiusbin);

    //Filtered tree
    TFile* fout = new TFile("Filtered_Pion_Output_trkr_nozero.root","RECREATE");
    TTree* filteree =(TTree*)intree->CopyTree("(Radiusbin[0]>0)&&(truth_Trkr)"); 
    
    //Selects entries for which first Radius bin is nonempty and event is wholly in trkr
    
    filteree->Write();
    fout->Close();
    delete fout;
}

Int_t* Radbinner(Double_t radarray[],Int_t arsize,Double_t binsize, Int_t radiusrange)
{
    Int_t radind = -99;
    Int_t* radbinarray = new Int_t[radiusrange];

    for(Int_t ii = 0; ii<arsize; ii++)
    {
        if(radarray[ii]>1e-12)
        {
            radind = floor(radarray[ii]/binsize);
            radbinarray[radind] += 1;
        }
    }

    return radbinarray;
}


Int_t FindRadMax(Int_t binarray[], Int_t arsize)
{
    Int_t radmax = -99;

    for(Int_t ii = 0; ii<(arsize-1); ii++)
    {
        if(binarray[ii] == 0 && binarray[ii+1] == 0){radmax = ii; break;};
    }

    if(radmax == -99){fprintf(stderr,"\n\n infinite \n\n");exit(1);};
    
    return radmax;

}


void RadAna()
{
    TFile* file = new TFile("Filtered_Pion_Output_trkr_nozero.root","UPDATE");
    TTree* tree = (TTree*)file->Get("tout");


    Int_t radiusrange;
    Int_t Radiusbin[200] = {0};
    Double_t Energybin[200] = {0.0};
    tree->SetBranchAddress("radiusrange",&radiusrange);
    tree->SetBranchAddress("Radiusbin",&Radiusbin);
    tree->SetBranchAddress("Energybin",&Energybin);

    Double_t rmax, enovflow;
    Int_t indmax, ovflow;
    Double_t Radiusbincut[200] = {0.0};
    Double_t Energybincut[200] = {0.0};

    TBranch *rmaxbranch = tree->Branch("rmax",&rmax,"rmax/D");
    TBranch *indmaxbranch = tree->Branch("indmax",&indmax,"indmax/I"); 
    TBranch *ovflowbranch = tree->Branch("ovflow",&ovflow,"ovflow/I");
    TBranch *enovflowbranch = tree->Branch("enovflow",&enovflow,"enovflow/D");


    TBranch *radcutbranch = tree->Branch("Radiusbincut",&Radiusbincut,"Radiusbincut[200]/D");
    TBranch *energycutbranch = tree->Branch("Energybincut",&Energybincut,"Energybincut[200]/D");

    TH2D* histo = new TH2D("histo","Correlation between number of clusters at rmax and ovflow",15,0,15,15,0,15);




    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        //for(Int_t ii = 0;ii<radiusrange;ii++){printf(" %d ",Radiusbin[ii]);};
        indmax = FindRadMax(Radiusbin,radiusrange);
        rmax = indmax*5.0; //in cm 
        //printf("%g \n",rmax);

        ovflow = 0;
        enovflow = 0;

        for(Int_t ii = 0; ii<200;ii++){Radiusbincut[ii] = 0;Energybincut[ii] = 0;};

        for(Int_t ii = 0; ii<indmax; ii++){Radiusbincut[ii] = Radiusbin[ii];Energybincut[ii] = Energybin[ii];};
        for(Int_t ii = indmax;ii<radiusrange;ii++){ovflow+=Radiusbin[ii];enovflow+=Energybin[ii];};/*if(Energybin[ii] > 0) printf("\n abcdefg %llu %d %g %g  \n",jentry,ii,enovflow,Energybin[ii]);};*/


//test 
//if(enovflow > 0) exit(1);


//Prints arrays to log 
        printf("\n");
        for(Int_t ii = 0; ii<indmax; ii++){printf(" %g ",Radiusbincut[ii]);};
        printf(" \t %d \t\t",ovflow);
        for(Int_t ii = 0; ii<indmax+6; ii++){printf(" %d ",Radiusbin[ii]);};
        printf("\t\t");
        for(Int_t ii = 0; ii<indmax; ii++){printf(" %g ", Energybincut[ii]);};


        histo->Fill(Radiusbincut[indmax-1],ovflow);

        rmaxbranch->Fill();
        indmaxbranch->Fill();
        radcutbranch->Fill();
        energycutbranch->Fill();
        ovflowbranch->Fill();
        enovflowbranch->Fill();
    }

    histo->SetXTitle("Number of clusters at rmax");
    histo->SetYTitle("Overflow Clusters");

    TH2D* normhisto = PlotHelp::NormalHist(histo);
    normhisto->SetStats(0);

    tree->Write("",TObject::kOverwrite);
    histo->Write("",TObject::kOverwrite);
    normhisto->Write("",TObject::kOverwrite);
    file->Save();
    file->Close();

    delete file;

}



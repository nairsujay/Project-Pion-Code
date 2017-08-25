#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"

#include <vector>

#include "PlotHelp.h"


void RadiusMean()
{
    TFile::Open("Filtered_Pion_Output_trkr_nozero.root");
    TTree* tree = (TTree*)gDirectory->Get("tout");
    TList* lout = new TList;

    Int_t secprocnum, arrsize;

    std::vector<double> *Energyarr = new std::vector<double>;
    std::vector<double> *Radiusarr = new std::vector<double>;

    Double_t rmax;
//Debug related
//------------------------------------------------
    //Double_t Energybin[200];
    //tree->SetBranchAddress("Energybin",&Energybin);
//--------------------------------------------------

    tree->SetBranchAddress("Energyunbin",&Energyarr);
    tree->SetBranchAddress("Radiusunbin",&Radiusarr);
    tree->SetBranchAddress("arrsize",&arrsize);
    tree->SetBranchAddress("rmax",&rmax);
    tree->SetBranchAddress("secondary_proc",&secprocnum);

    TH2D* trunchist60 = new TH2D("trunchist60","Truncated Mean Energy (60%) against process;Process;Truncated Mean Energy (MeV)",8,0,8,25,0,100);
    TH2D* trunchist90 = new TH2D("trunchist90","Truncated Mean Energy (90%) against process;Process;Truncated Mean Energy (MeV)",8,0,8,25,0,100);
    TH2D* radiushist = new TH2D("radiushist","Mean Radius against Process;Process;Mean Radius (cm)",8,0,8,20,0,100);
    TH2D* rwehist = new TH2D("rwehist","Weighted Mean Radius (by Energy) against Process;Process;Weighted Mean Radius (cm)",8,0,8,20,0,100);
    TH2D* ewrhist = new TH2D("ewrhist","Weighted Mean Energy (by Radius) against Process;Process;Weighted Mean Energy (MeV)",8,0,8,20,0,100);

    lout->Add(trunchist60);
    lout->Add(trunchist90);
    lout->Add(radiushist);
    lout->Add(rwehist);
    lout->Add(ewrhist);

    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        std::vector<double> Energy; //New vectors to contain only values within rmax
        std::vector<double> Radius;


        tree->GetEntry(jentry);
        for(Int_t ii = 0; ii<arrsize;ii++)
        {
            if(((*Radiusarr)[ii]<rmax))
            {
                Energy.push_back((*Energyarr)[ii]);
                Radius.push_back((*Radiusarr)[ii]);
            }
        }

        Int_t ensize = Energy.size();
        Int_t rsize = Radius.size();

        //if(ensize == 0)
        //{
            //printf("\n");
            //for(Int_t ii = 0; ii <arrsize;ii++){printf(" %g ",(*Radiusarr)[ii]);};
            //printf("\t\t %g ",rmax);
            //printf("\t\t %g ",Energybin[0]);
        //}

        if(ensize != rsize){fprintf(stderr,"\n\n \t FAIL \t \n\n");exit(1);};

        Int_t sortindex[ensize]; //Array for sorted indices

        TMath::Sort(ensize,Energy.data(),sortindex,kFALSE); //Generate sorted indices

        Int_t ind60 = TMath::Ceil(0.6*ensize);
        Int_t ind90 = TMath::Ceil(0.9*ensize);

        //Calculate truncated means;
        Double_t truncmean60 = 0.0;
        Double_t truncmean90 = 0.0;

        for(Int_t ii = 0; ii<ind60;ii++){truncmean60 += Energy[sortindex[ii]];};
        for(Int_t ii = 0; ii<ind90;ii++){truncmean90 += Energy[sortindex[ii]];};

        truncmean60 /= ind60;
        truncmean90 /= ind90;

        trunchist60->Fill(secprocnum,truncmean60);
        trunchist90->Fill(secprocnum,truncmean90);

        //Calculate Radial Mean
        Double_t radmean = 0.0;
        for(Int_t ii = 0; ii<ensize;ii++){radmean +=Radius[ii];};
        radmean /= ensize;

        radiushist->Fill(secprocnum,radmean);


        //Calculate Weighted means
        Double_t rwemean = 0.0;
        Double_t ewrmean = 0.0;
        Double_t energytot = 0.0;
        Double_t radiustot = 0.0;
        for(Int_t ii = 0; ii<ensize;ii++)
        {
            rwemean += Radius[ii]/Energy[ii];
            ewrmean += Energy[ii]/Radius[ii];
            energytot += Energy[ii];
            radiustot += Radius[ii];
        }

        rwemean = (rwemean*energytot)/ensize;
        ewrmean = (ewrmean*radiustot)/ensize;

        rwehist->Fill(secprocnum,rwemean);
        ewrhist->Fill(secprocnum,ewrmean);

        printf(" %d \t\t %d %d \t\t  %g %g \t\t %g %g %g \n",ensize,ind60,ind90,truncmean60,truncmean90,radmean,rwemean,ewrmean);
    }

   //loop to label histos 
    for(int jj=1;jj<=8;jj++)
    {
        rwehist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        ewrhist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        radiushist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        trunchist90->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        trunchist60->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

        rwehist->SetStats(0);
        ewrhist->SetStats(0);
        radiushist->SetStats(0);
        trunchist90->SetStats(0);
        trunchist60->SetStats(0);


    TH2D*trunchistnorm60 = PlotHelp::NormalHist(trunchist60);
    TH2D*trunchistnorm90 = PlotHelp::NormalHist(trunchist90);
    TH2D* radiushistnorm = PlotHelp::NormalHist(radiushist);
    TH2D* rwenorm        = PlotHelp::NormalHist(rwehist);
    TH2D* ewrnorm        = PlotHelp::NormalHist(ewrhist);

    lout->Add(trunchistnorm60);
    lout->Add(trunchistnorm90);
    lout->Add(radiushistnorm);
    lout->Add(ewrnorm);
    lout->Add(rwenorm);

    TFile* fout = new TFile("outplot/Mean.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
}


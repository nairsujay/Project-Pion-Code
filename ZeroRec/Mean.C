#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"

#include <vector>

#include "PlotHelp.h"


const Color_t NiceColorArray[8] = {kBlack,kBlue,kGreen,kYellow,kMagenta,kRed,kCyan,kGray};

const char *capproclabel[8] = {"Absorption","Multi-Production", "Charge Exchange","Double Charge Exchange","Inelastic", 
    "Rangeout","Decay","Conversion"};

void BasicMean()
{
    TFile::Open("Pion_Output_Recur.root");
    TTree* tree = (TTree*)gDirectory->Get("tout");
    TList* lout = new TList;

    Int_t secprocnum, arrsize;
    Int_t modsize;
    Bool_t truth_Trkr;

    std::vector<double> *Energyarr = new std::vector<double>;
    std::vector<double> *Radiusarr = new std::vector<double>;


    tree->SetBranchAddress("Energyunbin",&Energyarr);
    tree->SetBranchAddress("Radiusunbin",&Radiusarr);
    tree->SetBranchAddress("arrsize",&arrsize);
    tree->SetBranchAddress("modsize",&modsize);
    tree->SetBranchAddress("truth_Trkr",&truth_Trkr);
    tree->SetBranchAddress("secondary_proc",&secprocnum);

    TH2D* clusth = new TH2D("clusth","Number of Clusters by process;Process;Number of Clusters",8,0,8,60,0,60);
    TH2D* radh = new TH2D("radh","Radius Distribution by process;Process;Radius (cm)",8,0,8,20,0,100);
    TH2D* energh = new TH2D("energh","Energy Distribution by process;Process;Energy Deposited (MeV)",8,0,8,20,0,40);

    lout->Add(clusth);
    lout->Add(radh);
    lout->Add(energh);

    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(truth_Trkr) //Only analyse trkr events
        {
            clusth->Fill(secprocnum,arrsize);

            for(Int_t ii = 0;ii < arrsize;ii++)
            {
                energh->Fill(secprocnum,(*Energyarr)[ii],1.0/arrsize);
                radh->Fill(secprocnum,(*Radiusarr)[ii],1.0/arrsize);
            }
        }
    }

    //loop to label hists
    for(int jj=1;jj<=8;jj++)
    {
        clusth->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        radh->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        energh->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    clusth->SetStats(0);
    radh->SetStats(0);
    energh->SetStats(0);

    TH2D* clusthn = PlotHelp::NormalHist(clusth);
    TH2D* radhn = PlotHelp::NormalHist(radh);
    TH2D* energhn = PlotHelp::NormalHist(energh);


    lout->Add(clusthn);
    lout->Add(radhn);
    lout->Add(energhn);

    TFile* fout = new TFile("outplot/BasicMeans.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
}


void RadiusMean()
{
    TFile::Open("Pion_Output_Recur.root");
    TTree* tree = (TTree*)gDirectory->Get("tout");
    TList* lout = new TList;

    Int_t secprocnum, arrsize;
    Bool_t truth_Trkr;

    std::vector<double> *Energyarr = new std::vector<double>;
    std::vector<double> *Radiusarr = new std::vector<double>;

//Debug related
//------------------------------------------------
    //Double_t Energybin[200];
    //tree->SetBranchAddress("Energybin",&Energybin);
//--------------------------------------------------

    tree->SetBranchAddress("Energyunbin",&Energyarr);
    tree->SetBranchAddress("Radiusunbin",&Radiusarr);
    tree->SetBranchAddress("arrsize",&arrsize);
    tree->SetBranchAddress("secondary_proc",&secprocnum);
    tree->SetBranchAddress("truth_Trkr",&truth_Trkr);

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
        tree->GetEntry(jentry);
        if(truth_Trkr && arrsize > 0) //Only analyse trkr events
        {
            Int_t sortindex[arrsize]; //Array for sorted indices

            TMath::Sort(arrsize,Energyarr->data(),sortindex,kFALSE); //Generate sorted indices in ascending order

            Int_t ind60 = TMath::Ceil(0.6*arrsize);
            Int_t ind90 = TMath::Ceil(0.9*arrsize);

            //Calculate truncated means;
            Double_t truncmean60 = 0.0;
            Double_t truncmean90 = 0.0;

            for(Int_t ii = 0; ii<ind60;ii++){truncmean60 += (*Energyarr)[sortindex[ii]];};
            for(Int_t ii = 0; ii<ind90;ii++){truncmean90 += (*Energyarr)[sortindex[ii]];};

            truncmean60 /= ind60;
            truncmean90 /= ind90;

            trunchist60->Fill(secprocnum,truncmean60);
            trunchist90->Fill(secprocnum,truncmean90);

            //Calculate Radial Mean
            Double_t radmean = 0.0;
            for(Int_t ii = 0; ii<arrsize;ii++){radmean += (*Radiusarr)[ii];};
            radmean /= arrsize;

            radiushist->Fill(secprocnum,radmean);


            //Calculate Weighted means
            Double_t rwemean = 0.0;
            Double_t ewrmean = 0.0;
            Double_t energytot = 0.0;
            Double_t radiustot = 0.0;
            for(Int_t ii = 0; ii<arrsize;ii++)
            {
                rwemean += (*Radiusarr)[ii]/(*Energyarr)[ii];
                ewrmean += (*Energyarr)[ii]/(*Radiusarr)[ii];
                energytot += (*Energyarr)[ii];
                radiustot += (*Radiusarr)[ii];
            }

            rwemean = (rwemean*energytot)/arrsize;
            ewrmean = (ewrmean*radiustot)/arrsize;

            rwehist->Fill(secprocnum,rwemean);
            ewrhist->Fill(secprocnum,ewrmean);

            printf(" %d \t\t %d %d \t\t  %g %g \t\t %g %g %g \n",arrsize,ind60,ind90,truncmean60,truncmean90,radmean,rwemean,ewrmean);
        }
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


void ProjClus()
{
    TFile* fin = new TFile("outplot/BasicMeans.root","READ");

    TH2D* histo = (TH2D*)fin->Get("clusth");

    TH1D* phist[8];
    THStack* hs = new THStack("allhist","Distribution of Cluster Number");
    THStack* hs2 = new THStack("othhist","Distribution of Cluster Number");
    TLegend* leg = new TLegend(0.55,0.55,0.75,0.85);
    TLegend* leg2 = new TLegend(0.65,0.65,0.85,0.85);

    Double_t scale[8] = {0.0};

    TList* lout = new TList;

    TCanvas* c1 = new TCanvas("c1","",900,600);

    for(Int_t ii = 0; ii<7; ii++)
    {
        TH1D* tmph = (TH1D*)histo->ProjectionY(secproclabel[ii],ii+1,ii+1); 


        tmph->SetLineWidth(2);
        tmph->SetLineColor(NiceColorArray[ii]);
        scale[ii] = tmph->Integral(); 
        tmph->Scale(1/scale[ii]);

       phist[ii] = (TH1D*)tmph->Clone(); 

       hs->Add(phist[ii]);
       leg->AddEntry(phist[ii],capproclabel[ii]);

       delete tmph;
        
    }


    phist[7] = (TH1D*)phist[0]->Clone("othhist");
    phist[7]->Scale(scale[0]);

    for(Int_t ii = 1; ii<5;ii++)
    {
        phist[7]->Add(phist[ii],scale[ii]);
    }

    scale[7] = phist[7]->Integral();

    phist[7]->Scale(1.0/scale[7]);

    phist[7]->SetLineColor(kBlack);

    style::ResetStyle(leg);
    style::ResetStyle(leg2);


    hs2->Add(phist[5]);
    hs2->Add(phist[7]);
    leg2->AddEntry(phist[5],"Rangeout");
    leg2->AddEntry(phist[7],"Others");

    hs->Draw("NOSTACK");
    hs->GetXaxis()->SetTitle("Number of Clusters");
    hs->GetYaxis()->SetTitle("Frequency");
    leg->Draw();
    c1->Modified();
    c1->Print("outplot/meanplots/allprocs.png");

    hs2->Draw("NOSTACK");
    hs2->GetXaxis()->SetTitle("Number of Clusters");
    hs2->GetYaxis()->SetTitle("Frequency");
    leg2->Draw();
    c1->Modified();
    c1->Print("outplot/meanplots/othprocs.png");

    lout->Add(hs);
    lout->Add(hs2);

    TFile* fout = new TFile("outplot/meanplots/Clusters.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;

}


void CutPlot()
{
    TList* lout = new TList;

    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");

    TH1D* histo = new TH1D("histo","Processes before  and after cut (Area Normalized)",8,0,8);
    TH1D* hist = new TH1D("hist","Processes after cut (Area Normalized)",8,0,8);

    Int_t secprocnum, arrsize;

    tree->SetBranchAddress("arrsize",&arrsize);
    tree->SetBranchAddress("secondary_proc",&secprocnum);

    for(Long64_t jentry = 0;jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(secprocnum != 6)  ////Filter out decay events
        {
        histo->Fill(secprocnum);

        if(arrsize<4) hist->Fill(secprocnum);
        }
    }

    Double_t scal = hist->Integral();
    Double_t scalo = histo->Integral();


    hist->Scale(1/scal);
    histo->Scale(1/scalo);

    histo->SetLineColor(kBlue);
    hist->SetLineColor(kRed);

    hist->SetMaximum(0.6);
    histo->SetMaximum(0.6);

    hist->SetStats(0);
    histo->SetStats(0);

    TLegend* leg = new TLegend(0.25,0.65,0.45,0.90);
    leg->AddEntry(histo,Form("Default ( %g )", 1.0));
    leg->AddEntry(hist,Form("Cut ( %g )",scal/scalo));


    for(int jj=1;jj<=8;jj++)
    {
        hist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        histo->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    style::ResetStyle(leg);

    TCanvas* c1 = new TCanvas("c1","",900,600);
    histo->Draw();
    hist->Draw("SAME");
    leg->Draw();
    c1->Print("outplot/cutplot.png");

    lout->Add(hist);
    lout->Add(histo);
    TFile* fout = new TFile("outplot/cutplot.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete fout;
    delete c1;
}


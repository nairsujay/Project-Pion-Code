#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

#include <vector>

#include "PlotHelp.h"



void CorrPlot(const TString varx,const TString vary, const TString Title,
              Int_t nbinsx,Double_t xmin,Double_t xmax, 
              Int_t nbinsy,Double_t ymin,Double_t ymax,
              const TString outdir, Bool_t kWeight)
{
    TFile* fin = new TFile("Pion_Output_Recur.root","READ");
    TTree* tree = (TTree*)fin->Get("tout");

    std::vector<Double_t> *Xarray = new std::vector<Double_t>;
    std::vector<Double_t> *Yarray = new std::vector<Double_t>;

    TH2D* lohist[8];
    TH2D* hihist[8];

    TList* lout = new TList;

    TCanvas* c1 = new TCanvas("c1","",900,600);

    for(Int_t ii = 0; ii < 8; ii++)
    {
        lohist[ii] = new TH2D(Form("lo%s",secproclabel[ii]),Title,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
        hihist[ii] = new TH2D(Form("hi%s",secproclabel[ii]),Title,nbinsx,xmin,xmax,nbinsy,ymin,ymax);

        lohist[ii]->SetMarkerColor(kBlue);
        hihist[ii]->SetMarkerColor(kRed);

        lohist[ii]->SetLineColor(kBlue);
        hihist[ii]->SetLineColor(kRed);

        lohist[ii]->SetStats(0);
        hihist[ii]->SetStats(0);

        lohist[ii]->GetXaxis()->SetTitleOffset(1.8);
        lohist[ii]->GetYaxis()->SetTitleOffset(1.8);


        hihist[ii]->GetXaxis()->SetTitleOffset(1.8);
        hihist[ii]->GetYaxis()->SetTitleOffset(1.8);
    }

    Int_t secprocnum, arrsize;
    Bool_t truth_Trkr;

    tree->SetBranchAddress(varx,&Xarray);
    tree->SetBranchAddress(vary,&Yarray);
    tree->SetBranchAddress("secondary_proc",&secprocnum);
    tree->SetBranchAddress("arrsize",&arrsize);
    tree->SetBranchAddress("truth_Trkr",&truth_Trkr);



    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(truth_Trkr)
        {
            for(Int_t ii = 0; ii<arrsize;ii++)
            {

                if(kWeight)
                {
                    if(arrsize<4) lohist[secprocnum]->Fill((*Xarray)[ii],(*Yarray)[ii]);
                    else hihist[secprocnum]->Fill((*Xarray)[ii],(*Yarray)[ii]);
                }

                if(!kWeight)
                {
                    if(arrsize<4) lohist[secprocnum]->Fill((*Xarray)[ii],(*Yarray)[ii],1.0/arrsize);
                    else hihist[secprocnum]->Fill((*Xarray)[ii],(*Yarray)[ii],1.0/arrsize);
                }
            }       
        }
    }

    for(Int_t ii = 0; ii < 8; ii++)
    {        
        lout->Add(lohist[ii]);
        lout->Add(hihist[ii]);


        lohist[ii]->Draw("CONT");
        c1->Print(Form("outplot/%s/lo%s.png",outdir.Data(),secproclabel[ii]));
        hihist[ii]->Draw("CONT");
        c1->Print(Form("outplot/%s/hi%s.png",outdir.Data(),secproclabel[ii]));


        lohist[ii]->Draw("BOX");
        hihist[ii]->Draw("BOX SAME");
        c1->Print(Form("outplot/%s/all%s.png",outdir.Data(),secproclabel[ii]));
    }

    TFile* fout = new TFile(Form("outplot/%s/Corrs.root",outdir.Data()),"RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
    delete c1;
}

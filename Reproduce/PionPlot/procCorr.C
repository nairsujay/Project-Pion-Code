#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"


void procCorr()
{

 //   Int_t proc_num = -99;
 //   Int_t secproc_num = -99;
   // TFile::Open("Pion_Output_MINERvA.root");
    TTree* tree = (TTree*)gDirectory->Get("tout");
    //tree->SetBranchAddress("proc",&proc_num);
    //tree->SetBranchAddress("secondary_proc",&secproc_num);

    TH1I* prochist = new TH1I("prochist","Process Frequencies",12,0,12);
    TH1I* secprochist = new TH1I("secprochist","Secondary-Process Frequencies",8,0,8);
    TH2I* proccorrhist = new TH2I("proccorrhist","Process and Secondary Process Correlation"
                                  ,12,0,12,8,0,8);
    const char *proclabel[12] = {"hadElastic","NeutronInelastic","ProtonInelastic",
                                 "PionMinusInelastic","PionPlusInelastic",
                                 "CHIPS","Decay","conv", "compt", 
                                 "Primary","Other","Rangeout"};
    const char *secproclabel[8] = {"abs","multiprod", "cex","dbl_cex","inel", 
                                    "rangeout","decay","conv"};

    //main loop for filling histogram -----Made Unnecessary
    /*for(int ii=0;ii<tree->GetEntries();ii++)*/
    //{
       //tree->GetEntry(ii);
       ////cout<<proc_num<<"\t"<<secproc_num<<endl;
       ////prochist->Fill(proc_num);
       ////secprochist->Fill(secproc_num);
       ////proccorrhist->Fill(proc_num,secproc_num);

    /*}*/

    

    //TESTING
    TCanvas* c2 = new TCanvas("c2","");
    tree->Draw("proc>>+prochist");
    tree->Draw("secondary_proc>>+secprochist");
    tree->Draw("secondary_proc:proc>>+proccorrhist");
  

//Styling
    prochist->SetStats(0);
    secprochist->SetStats(0);
    proccorrhist->SetStats(0);
    
    prochist->SetFillColor(4);
    prochist->SetXTitle("Process");
    prochist->SetYTitle("Events");

    secprochist->SetFillColor(4);
    secprochist->SetXTitle("Secondary Process");
    secprochist->SetYTitle("Events");

    proccorrhist->SetXTitle("Process");
    proccorrhist->SetYTitle("Secondary Process");



    //Loops to label bins
    for(int jj=1;jj<=12;jj++)
    {
        prochist->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
        proccorrhist->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
    }

    for(int jj=1;jj<=8;jj++)
    {
        secprochist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        proccorrhist->GetYaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    TFile* fout = new TFile("outplot/proccorr/proccorr.root","RECREATE");
    prochist->Write();
    secprochist->Write();
    proccorrhist->Write();
    fout->Close();
    TCanvas* c1= new TCanvas("c1","proc");
    c1->SetGrid();
    prochist->Draw();
    c1->Print("outplot/proccorr/proc.png");
    secprochist->Draw();
    c1->Print("outplot/proccorr/secproc.png");
    proccorrhist->Draw("COLZ");
    c1->Print("outplot/proccorr/proccorr.png");
}


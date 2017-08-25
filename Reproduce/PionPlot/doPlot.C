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


//const TString gcut("(ecal10==0)"); //gets rid of """false""" zeroes BE CAREFUL WITH THIS
const TString gcut("(trkr10!=0)&&(ecal == 0)"); //gets rid of """false""" zeroes [[[STILL IN TESTING]]]


//Reproduces plots in paper
void ReproducePlot()
{

    TList* lout = new TList ;

    const TString ecalcut("&&trkr==0"); // only trkr tracks
    const TString trkrcut("&&ecal==0"); // only ecal tracks
    const TString trkrecal(""); //let's all in ?
    const TString bothonly("&&((trkr>0)&&(ecal>0))");

    const TString Proc[8] = {"Absorption","Multi-Production","Charge Exchange","Double Charge Exchange","Inelastic","Rangeout","Decay","all"};
    const TString Cut[8] = {"secondary_proc==0","secondary_proc==1","secondary_proc==2",
                            "secondary_proc==3","secondary_proc==4","secondary_proc==5",
                            "secondary_proc==6","pdg==211"};
    // pdg cut does not remove any data, 
    //it's purely for syntax as Plot1DModule requires a cut

    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");

    TCanvas* c1 = new TCanvas("c1","",600,400);

    for(Int_t ii=0;ii<8;ii++) //Nested loops iterate over process and radius
    {
        for(Int_t jj=10;jj<=30;jj+=5)
        {        
            PlotHelp::NicePlot(tree,c1,Proc[ii],Form("trkr%d",jj),Cut[ii]+trkrcut,60,0,160,lout,"outplot/reproduce",kFALSE,jj);

            PlotHelp::NicePlot(tree,c1,Proc[ii],Form("ecal%d",jj),Cut[ii]+ecalcut,60,0,160,lout,"outplot/reproduce",kFALSE,jj);
        }
    }

    delete c1;


    TFile* fout =new TFile("outplot/reproduce/EnergyDepositbyProcess.root","RECREATE");
    lout->Write();
    fout->Close();
    delete lout;
    delete fout;
}

//3 different cuts at trkr30 
void CutPlot(Double_t thresh_1,Double_t thresh_2,Double_t thresh_3) 
{
    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");

    TList* lout = new TList;
    Double_t threshold[4] = {0,thresh_1,thresh_2,thresh_3};
    Double_t procnorm[4] = {-1,-1,-1,-1}; //Contains number of entries for each cut
    TH1D* prochist[4];
    TH1D* secprochist[4];
    TH2D* proccorrhist[4];
    TString allcut[4]; //TString array will contain cuts
    allcut[0] = TString("");

    for(Int_t ii = 1; ii<4;ii++){allcut[ii] = TString(Form("(trkr30<%g)&&%s",threshold[ii],gcut.Data()));}

    //loops through cuts and plots them individually
    for(Int_t ii = 0; ii<4; ii++)
    {
        TCanvas* c1 = new TCanvas("c1","",1280,900);
        c1->SetGrid();
        procnorm[ii] = PlotHelp::ProcPlot(tree,c1,allcut[ii],
                "cutproccorr",prochist[ii],secprochist[ii],proccorrhist[ii],ii);

        lout->Add(prochist[ii]);
        lout->Add(secprochist[ii]);
        lout->Add(proccorrhist[ii]);

        delete c1;
    }

    //Adds all hists to a THStack to be plotted together
    TCanvas* c1 = new TCanvas("c1","",1280,900);
    PlotHelp::StackProc(prochist,c1,lout,allcut,"cutproccorr",procnorm,4,kFALSE);
    PlotHelp::StackProc(secprochist,c1,lout,allcut,"cutproccorr",procnorm,4,kTRUE);
    delete c1;


    TFile* fout = new TFile("outplot/cutproccorr/CutProc.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();
    
    delete fout;
    delete lout;
}




//Cuts across multiple radii or just one...
void MultiCut(const TString cut_1,const TString cut_2)
{
    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");

    TList* lout = new TList;
    Double_t procnorm[3] = {-1,-1,-1}; //array of entries of histograms
    TH1D* prochist[3];
    TH1D* secprochist[3];
    TH2D* proccorrhist[3];
    TString allcut[3] = {gcut,cut_1+"&&"+gcut,cut_2+"&&"+gcut}; //Array of cuts

    TCanvas* c1 = new TCanvas("c1","",1280,900);

    //loops thorugh cuts
    for(Int_t ii = 0; ii<3 ; ii++)
    {
        procnorm[ii] = PlotHelp::ProcPlot(tree,c1,allcut[ii],"multicut",prochist[ii],
                secprochist[ii],proccorrhist[ii],ii);
        //printf("\n\n%g\t%d\t%d\n\n",procnorm[ii],prochist[ii]->GetSize(),secprochist[ii]->GetSize());
        lout->Add(prochist[ii]);
        lout->Add(secprochist[ii]);
        lout->Add(proccorrhist[ii]);
    }   

    //Stacks histograms
    PlotHelp::StackProc(prochist,c1,lout,allcut,"multicut",procnorm,3,kFALSE);
    PlotHelp::StackProc(secprochist,c1,lout,allcut,"multicut",procnorm,3,kTRUE);

    TFile* fout = new TFile("outplot/multicut/MultiCutProc.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete fout;
    delete lout;
}





//Creates 2D histogram of energy deposited at different radii
void TrkrCorr(Int_t proc)
{

    TFile* fin = new TFile("Pion_Output_MINERvA.root", "READ");
    TTree* tree =(TTree*)fin->Get("tout");

    TList* lout = new TList;
    TH2D* hout[4];

//loops through different radii
    for(Int_t ii=10;ii<30;ii += 5)
    {
        TCanvas* c1 =new TCanvas("c1","");
        TH2D* tmph = new TH2D("tmph","hist",100,0,300,100,0,300);
        tree->Draw(Form("trkr%d:trkr%d>>+tmph",ii+5,ii),gcut+Form("&&(secondary_proc==%d)",proc));
        TH2D* tmpnormh = PlotHelp::NormalHist(tmph); //Normalizes by slice

        hout[(ii/5)-2] =(TH2D*)tmpnormh->Clone(Form("trkr%d%d_%d",ii,ii+5,proc));
        style::ResetStyle(hout[(ii/5)-2]);
        lout->Add(hout[(ii/5)-2]);

        delete tmpnormh;
        delete tmph;
        delete c1;
    }

    TFile* fout = new TFile("outplot/Trkrcorr.root","UPDATE");
    lout->Write();
    fout->Save();
    fout->Close();
}

   

void ZeroPlot()
{
    const TString zeroecal("&&(!((ecal10 == 0)&&(ecal15 == 0)&&(ecal20 == 0)&&(ecal25 == 0)&&(ecal30 == 0)))");
    const TString zerotrkr("&&(!((trkr10 == 0)&&(trkr15 == 0)&&(trkr20 == 0)&&(trkr25 == 0)&&(trkr30 == 0)))");
    const TString nozeroecal("&&(trkr == 0)&&(ecal10 !=0)&&(ecal15 != 0)&&(ecal20 != 0)&&(ecal25 != 0)&&(ecal30 != 0)");
    const TString nozerotrkr("&&(ecal == 0)&&(trkr10 !=0)&&(trkr15 != 0)&&(trkr20 != 0)&&(trkr25 != 0)&&(trkr30 != 0)");
    const TString ecalcut("&&(trkr==0)&&(ecal10 != 0)");
    const TString trkrcut("&&(ecal==0)&&(trkr10 != 0)");
    const TString origtrkrcut("&&(ecal==0)&&(trkr10 + trkr15 + trkr20 + trkr25 + trkr30 + trkr !=0)");
    const TString origecalcut("&&(trkr==0)&&(ecal10 + ecal15 + ecal20 + ecal25 + ecal30 + ecal !=0)");
    


    TList* lout = new TList ;
    const TString Proc[8] = {"Absorption","multiprod","cex","dbl_cex","Inelastic","Rangeout","decay","all"};
    const TString Cut[8] = {"secondary_proc==0","secondary_proc==1","secondary_proc==2",
                            "secondary_proc==3","secondary_proc==4","secondary_proc==5",
                            "secondary_proc==6","pdg==211"};
    const TString sProc[8] = {"sabs","smultiprod","scex","sdbl_cex","sinel","srangeout","sdecay","sall"};
    
    // pdg cut does not remove any data, 
    //it's purely for syntax as Plot1DModule requires a cut

    TFile* fin = new TFile("Pion_Output_MINERvA.root", "READ");
    TTree* tree =(TTree*)fin->Get("tout");


    TCanvas* c1 = new TCanvas("c1","",600,400);


    for(Int_t ii=0;ii<8;ii++) //Nested loops iterate over process and radius
    {
        for(Int_t jj=10;jj<=30;jj+=5)
        {
            PlotHelp::NicePlot(tree,c1,Proc[ii],Form("trkr%d",jj),Cut[ii]+trkrcut,60,0,160,lout,"outplot/zerocut",kTRUE,jj);

            PlotHelp::NicePlot(tree,c1,Proc[ii],Form("ecal%d",jj),Cut[ii]+ecalcut,60,0,160,lout,"outplot/zerocut",kTRUE,jj);

            //PlotHelp::Plot1DModule(tree,c1,sProc[ii],Form("trkr%d",jj),Cut[ii]+nozerotrkr,60,0,160,lout,"outplot/strictzerocut");

            //PlotHelp::Plot1DModule(tree,c1,sProc[ii],Form("ecal%d",jj),Cut[ii]+nozeroecal,60,0,160,lout,"outplot/strictzerocut");

        }
    }

    delete c1;


    TFile* fout =new TFile("outplot/zerocut/EnergyDepositbyProcess.root","RECREATE");
    lout->Write();


    fout->Close();
    delete lout;
    delete fout;
}



void TrkrRadCorr()
{
    const TString Proc[8] = {"abs","multiprod","cex","dbl_cex","inel","rangeout","decay","conv"};

    Int_t secproc = -99;
    Double_t trkr10, trkr15, trkr20, trkr25, trkr30;
    Double_t ecal=0;

    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");
    tree->SetBranchAddress("trkr10",&trkr10);
    tree->SetBranchAddress("trkr15",&trkr15);
    tree->SetBranchAddress("trkr20",&trkr20);
    tree->SetBranchAddress("trkr25",&trkr25);
    tree->SetBranchAddress("trkr30",&trkr30);
    tree->SetBranchAddress("ecal",&ecal);
    tree->SetBranchAddress("secondary_proc",&secproc);

    TH2D* secproctrkrhist[8];

    TList* lout = new TList;
    TCanvas* c1 = new TCanvas("c1","",900,600);

    //Main Loop over Processes
    for(Int_t ii = 0; ii<8; ii++)
    {
        TH2D* tmphist = new TH2D("tmph",Proc[ii],5,10,35,60,0,160);
        for(Int_t jj = 0; jj<tree->GetEntries();jj++)
        {
            tree->GetEntry(jj);
            if((secproc == ii)&&(trkr10 != 0)&&(ecal == 0))
            {
                tmphist->Fill(10,trkr10);
                tmphist->Fill(15,trkr15);
                tmphist->Fill(20,trkr20);
                tmphist->Fill(25,trkr25);
                tmphist->Fill(30,trkr30);
            }
        }

        tmphist->SetXTitle("Radius from Track Endpoint (cm)");
        tmphist->SetYTitle("cumulative Energy Deposited (MeV)");
        TH2D* tmpnormhist = PlotHelp::NormalHist(tmphist);
        style::ResetStyle(tmpnormhist);
        tmpnormhist->GetXaxis()->SetTitleOffset(1.0);
        tmpnormhist->GetYaxis()->SetTitleOffset(1.0);
        tmpnormhist->SetStats(0);
        
        tmpnormhist->Draw("COLZ");
        c1->Print(Form("outplot/trkrcorr/%s.png",Proc[ii].Data()));
        secproctrkrhist[ii] =(TH2D*)tmpnormhist->Clone(Proc[ii]);
        lout->Add(secproctrkrhist[ii]);

        delete tmphist;
    }
TFile* fout = new TFile("outplot/trkrcorr/proctrkr.root","RECREATE");
lout->Write();
fout->Save();
fout->Close();

delete fout;
delete lout;
delete c1;
}


void EcalRadCorr()
{
    const TString Proc[8] = {"abs","multiprod","cex","dbl_cex","inel","rangeout","decay","conv"};

    Int_t secproc = -99;
    Double_t ecal10, ecal15, ecal20, ecal25, ecal30;
    Double_t trkr=0;

    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");
    tree->SetBranchAddress("ecal10",&ecal10);
    tree->SetBranchAddress("ecal15",&ecal15);
    tree->SetBranchAddress("ecal20",&ecal20);
    tree->SetBranchAddress("ecal25",&ecal25);
    tree->SetBranchAddress("ecal30",&ecal30);
    tree->SetBranchAddress("trkr",&trkr);
    tree->SetBranchAddress("secondary_proc",&secproc);

    TH2D* secprocecalhist[8];

    TList* lout = new TList;
    TCanvas* c1 = new TCanvas("c1","",900,600);

    //Main Loop over Processes
    for(Int_t ii = 0; ii<8; ii++)
    {
        TH2D* tmphist = new TH2D("tmph",Proc[ii],5,10,35,60,0,160);
        for(Int_t jj = 0; jj<tree->GetEntries();jj++)
        {
            tree->GetEntry(jj);
            if((secproc == ii)&&(ecal10 != 0)&&(trkr == 0))
            {
                tmphist->Fill(10,ecal10);
                tmphist->Fill(15,ecal15);
                tmphist->Fill(20,ecal20);
                tmphist->Fill(25,ecal25);
                tmphist->Fill(30,ecal30);
            }
        }

        tmphist->SetXTitle("Radius (cm)");
        tmphist->SetYTitle("Energy (MeV)");
        TH2D* tmpnormhist = PlotHelp::NormalHist(tmphist);
        style::ResetStyle(tmpnormhist);
        tmpnormhist->GetXaxis()->SetTitleOffset(0.8);
        tmpnormhist->GetYaxis()->SetTitleOffset(0.9);
        tmpnormhist->SetStats(0);
        
               tmpnormhist->Draw("COLZ");
        c1->Print(Form("outplot/ecalcorr/%s.png",Proc[ii].Data()));
        secprocecalhist[ii] =(TH2D*)tmpnormhist->Clone(Proc[ii]);
        lout->Add(secprocecalhist[ii]);

        delete tmphist;
    }
TFile* fout = new TFile("outplot/ecalcorr/procecal.root","RECREATE");
lout->Write();
fout->Save();
fout->Close();

delete fout;
delete lout;
delete c1;
}

void ZeroProc()
{
const TString invgcut("(ecal==0)&&(trkr10==0)");

    TList* lout = new TList;

    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");

    TH1D* zeroproc = 0x0;
    TH1D* zerosec = 0x0;
    TH2D* zerocorr = 0x0;

    TCanvas* c1 = new TCanvas("c1","",1290,900);

    PlotHelp::ProcPlot(tree,c1,invgcut,"invcut",zeroproc,zerosec,zerocorr,1);

    zerosec->SetTitle("Cut Entries by Secondary Process (Area Normalized);Process;Frequency");
    zerosec->GetXaxis()->SetTitleOffset(0.8);
    zerosec->GetYaxis()->SetTitleOffset(0.8);

    lout->Add(zeroproc);
    lout->Add(zerosec);
    lout->Add(zerocorr);
    
    TFile* fout = new TFile("outplot/invcut/ZeroCutProc.root","RECREATE");
    lout->Write();
    fout->Close();

    delete lout;
    delete fout;
}


void TrkrvsEcal()
{

    TList* lout = new TList ;

    const TString ecalcut("&&(trkr + trkr10 + trkr15 + trkr20 + trkr25 + trkr30 == 0) "); // only trkr tracks
    const TString trkrcut("&&(ecal + ecal10 + ecal15 + ecal20 + ecal25 + ecal30 == 0)"); // only ecal tracks
    const TString trkrecal(""); //let's all in ?
    const TString bothonly("&&((trkr>0)&&(ecal>0))");

    const TString Proc[8] = {"Absorption","multiprod","cex","dbl_cex","Inelastic","Rangeout","decay","all"};
    const TString Cut[8] = {"secondary_proc==0","secondary_proc==1","secondary_proc==2",
                            "secondary_proc==3","secondary_proc==4","secondary_proc==5",
                            "secondary_proc==6","pdg==211"};
    const TString names[4] = {"ecalcut","trkrcut","trkrecal","bothonly"};
    const TString znames[4] = {"zecalcut","ztrkrcut","ztrkrecal","zbothonly"};
    
    const TString encut[4] = {ecalcut,trkrcut,trkrecal,bothonly};
    const TString zencut[4] = {ecalcut+"&&(ecal10>0)",trkrcut+"&&(trkr10>0)",trkrecal+"&&(ecal10+trkr10>0)",bothonly+"&&(ecal10+trkr10>0)"};

    // pdg cut does not remove any data, 
    //it's purely for syntax as Plot1DModule requires a cut

    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");

    TCanvas* c1 = new TCanvas("c1","",600,400);
    for(Int_t ind1 = 0; ind1 < 4; ind1++)
    {
        for(Int_t ii=0;ii<8;ii++) //Nested loops iterate over process and radius
        {
            for(Int_t jj=10;jj<=30;jj+=5)
            {        
                PlotHelp::Plot1DModule(tree,c1,Proc[ii]+"_"+names[ind1],Form("trkr%d",jj),Cut[ii]+encut[ind1],60,0,160,lout,"outplot/trkrvsecal");

                PlotHelp::Plot1DModule(tree,c1,Proc[ii]+"_"+names[ind1],Form("ecal%d",jj),Cut[ii]+encut[ind1],60,0,160,lout,"outplot/trkrvsecal");

                PlotHelp::Plot1DModule(tree,c1,Proc[ii]+"_"+znames[ind1],Form("trkr%d",jj),Cut[ii]+zencut[ind1],60,0,160,lout,"outplot/trkrvsecal");

                PlotHelp::Plot1DModule(tree,c1,Proc[ii]+"_"+znames[ind1],Form("ecal%d",jj),Cut[ii]+zencut[ind1],60,0,160,lout,"outplot/trkrvsecal");

            }
        }
    }


    TFile* fout = new TFile("outplot/trkrvsecal/AATrkrvsEcal.root","RECREATE");
    lout->Write();
    fout->Close();


    delete c1;
    delete fout;
    delete lout;
}

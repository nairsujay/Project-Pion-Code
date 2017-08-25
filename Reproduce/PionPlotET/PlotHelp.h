#ifndef _PLOTHELP_H_
#define _PLOTHELP_H_


#include "style.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TList.h"
#include "THStack.h"
#include "stdlib.h"

//Global Variables
const Double_t EPSILON = 1e-12;
const char *proclabel[12] = {"hadElastic","NeutronInelastic","ProtonInelastic",
    "PionMinusInelastic","PionPlusInelastic",
    "CHIPS","Decay","conv", "compt", 
    "Primary","Other","Rangeout"};
const char *secproclabel[8] = {"abs","multiprod", "cex","dbl_cex","inel", 
    "rangeout","decay","conv"};

const Color_t ColorArray[4] = {kBlack,kBlue,kRed,kGreen};

//Class to help with plots 

class PlotHelp
{
    public:

        //Helpful Methods

        //Lifted from NeutrinoTools.h with minor stylistic changes
        static  void Plot1DModule(TTree *tt, TCanvas *c1, const TString tag, 
                     const TString var,const TString allcut, const Int_t nbins,
                     const Double_t xmin, const Double_t xmax, 
                     TList *list, const TString outdir,Int_t jj);

        //Lifted from NeutrinoTools.h
        static  TH2D* NormalHist(TH2D* hraw);

        //Draws a 1Dd histogram of secondary & primary processes & a 2D histogram of correlation and labels bins
        static  Double_t ProcPlot(TTree* tree,TCanvas* c1, const TString allcut,
                         const TString outdir,TH1D* &proc, TH1D* &secproc, TH2D* &proccorr, 
                         Int_t Color);

        //Creates a "name" string (histogram or file) from a given cut
        static  TString Cut2Name(const TString allcut);

        //Creates a THStack from an array of histograms and a legend
        static  void StackProc(TH1D* histo[],TCanvas* c1,TList* lout, 
                     const TString allcut[],const TString outdir,Double_t nor[],
                     Int_t nhists,Bool_t kSec);

        //File opener
        static  TTree* InputFile(const TString fname,const TString tname);

        //Nice plotter
        static  void NicePlot(TTree *tt, TCanvas *c1, const TString title, const TString var, 
                     const TString allcut, const Int_t nbins, const Double_t xmin, 
                     const Double_t xmax, TList *list, const TString outdir, 
                     Bool_t kZERO, Int_t ii);

};

TTree* PlotHelp::InputFile(const TString fname,const TString tname)
{
    TFile* fin = new TFile(fname,"READ");
    TTree* tree =(TTree*)fin->Get(tname);

    return tree;
}


void PlotHelp::NicePlot(TTree *tt, TCanvas *c1, const TString title, const TString var, 
        const TString allcut, const Int_t nbins, const Double_t xmin, 
        const Double_t xmax, TList *list, const TString outdir, 
        Bool_t kZERO, Int_t ii)
{
    gPad->SetGrid();
    gStyle->SetHistMinimumZero();

    TString hname=title+var;

    TH1D *hh=new TH1D(hname,title,nbins, xmin, xmax); 
    if(list){
        list->Add(hh);
    }

    tt->Draw(var+">>"+hname, allcut);
    printf("Plot1DModule title %s var %s allcut %s\n", title.Data(), var.Data(), allcut.Data());

    hh->SetMaximum(hh->GetBinContent(hh->GetMaximumBin())*1.05);

    //if(kZERO) hh->SetMaximum(8e3);
    //if(!kZERO) hh->SetMaximum(6e4);

    hh->SetStats(0);


    if(var.Contains("trkr"))
    {
        if(kZERO) hh->SetTitle(Form(" %s  (Zero Suppressed)",title.Data()));
        
        if(!kZERO) hh->SetTitle(Form(" %s ",title.Data()));
        
    }

    if(var.Contains("ecal"))
    {
        hh->SetTitle(title+Form(" - Energy Deposited at Ecal Radius = %d",ii));
    }

    hh->SetXTitle("Cumulative Energy Deposited (MeV)");
    hh->SetYTitle("Events");
    hh->GetYaxis()->SetTitleOffset(1.4);

    
    hh->SetFillColor(kRed);
    hh->SetLineColor(kRed); hh->SetLineWidth(2); hh->Draw();
    c1->Print(Form("%s/%s_%s.png",outdir.Data(), title.Data(),var.Data()));

    if(!list){
        delete hh;
    }
}



void PlotHelp::Plot1DModule(TTree *tt, TCanvas *c1, const TString tag, const TString var, 
        const TString allcut, const Int_t nbins, const Double_t xmin, 
        const Double_t xmax, TList *list, const TString outdir,Int_t jj)
      
{
    gPad->SetGrid();
    gStyle->SetHistMinimumZero();

    TString hname=tag;

    const TString Title[8] = {"Absorption","Multi-Production","Charge Exchange","Double Charge Exchange","Inelastic","Rangeout","Decay","all"};

    const Double_t maxh=1.1;
    TH1D *hh=new TH1D(hname,Title[jj],nbins, xmin, xmax); 
    if(list){
        list->Add(hh);
    }

    tt->Draw(var+">>"+hname, allcut);
    printf("Plot1DModule tag %s var %s allcut %s\n", tag.Data(), var.Data(), allcut.Data());

    hh->SetMaximum(hh->GetBinContent(hh->GetMaximumBin())*maxh);
    hh->SetStats(0); 
    hh->SetFillColor(kRed);
    hh->SetLineColor(kRed); hh->SetLineWidth(2); hh->Draw();
    c1->Print(Form("%s/%s.png",outdir.Data(), tag.Data()));

    if(!list){
        delete hh;
    }
}

//Normalizes a 2D hist (Can be overloaded to TH2I etc.) 
TH2D* PlotHelp::NormalHist(TH2D* hraw)
{
    TH2D* hh = (TH2D*)hraw->Clone(Form("%snor",hraw->GetName())); //Makes clone of input histo
    hh->Scale(0); //Clears clone

    //Finds bounds
    const Int_t x0 = hh->GetXaxis()->GetFirst();
    const Int_t x1 = hh->GetXaxis()->GetLast();
    const Int_t y0 = hh->GetYaxis()->GetFirst();
    const Int_t y1 = hh->GetYaxis()->GetLast();

    Double_t hmax = -1e10;
    Double_t hmin = 1e10;
    Double_t nent = 0;

    //Loops through slices and projects to Y 
    for(Int_t ix=x0; ix<=x1; ix++)
    {
        TH1D * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d",hh->GetName(), ix),ix,ix);
        const Double_t tot = sliceh->GetEntries();

        TH1D * pdfh=0x0;

        //Gets bin content of each bin along slice and finds max bin
        if(tot>EPSILON)
        {
            nent += tot;
            Double_t imax = -999;
            imax = sliceh->GetBinContent(sliceh->GetMaximumBin());

            //Rescales slice according to max bin content
            for(Int_t iy=y0; iy<=y1; iy++)
            {
                const Double_t cont = sliceh->GetBinContent(iy)/imax;
                if(tot>EPSILON && cont>0)
                {
                    hh->SetBinContent(ix,iy,cont);
                    if(cont>hmax) hmax = cont;
                    if(cont<hmin) hmin = cont;
                }
            }
        }

        delete pdfh;
        delete sliceh;
    }

    hh->SetEntries(nent);
    hh->SetMinimum(0.99*hmin);
    hh->SetMaximum(1.1*hmax);

    return hh;
}


//Plots processess and labels them etc. Returns number of entries for given cut
Double_t PlotHelp::ProcPlot(TTree* tree,TCanvas* c1, const TString allcut,const TString outdir
        , TH1D* &proc, TH1D* &secproc, TH2D* &proccorr, Int_t Color) 
{
    Double_t procnorm = -99;
    Double_t checknorm = -99;


    TH1D* prochisttmp = new TH1D("prochisttmp","Process Frequencies",12,0,12);
    TH1D* secprochisttmp = new TH1D("secprochisttmp","Secondary-Process Frequencies",8,0,8);
    TH2D* proccorrhisttmp = new TH2D("proccorrhisttmp",
            "Process and Secondary Process Correlation",
            12,0,12,8,0,8);


    tree->Draw("proc>>+prochisttmp",allcut);
    tree->Draw("secondary_proc>>+secprochisttmp",allcut);
    tree->Draw("secondary_proc:proc>>+proccorrhisttmp",allcut);


    //Styling

    style::ResetStyle(prochisttmp);
    style::ResetStyle(secprochisttmp);
    style::ResetStyle(proccorrhisttmp);

    prochisttmp->SetStats(0);
    secprochisttmp->SetStats(0);
    proccorrhisttmp->SetStats(0);

    prochisttmp->GetXaxis()->SetTickLength(0);
    secprochisttmp->GetXaxis()->SetTickLength(0);

    prochisttmp->SetLineWidth(2);
    secprochisttmp->SetLineWidth(2);

    prochisttmp->SetLineColor(ColorArray[Color]);
    secprochisttmp->SetLineColor(ColorArray[Color]);


    procnorm = prochisttmp->GetEntries();
    checknorm = secprochisttmp->GetEntries();

    printf(Form("\n\n\t Rangeout = %g \t\n\n",secprochisttmp->GetBinContent(6)));


    //Exit if secproc and proc have different entries--should be unlikely/impossible
    if(procnorm != checknorm){printf("\n\ndifferent norm\n\n");exit(1);};


    //normalize 1D histos by area and style them further
    prochisttmp->Scale(1.0/procnorm);
    secprochisttmp->Scale(1.0/procnorm);

    prochisttmp->SetMaximum(prochisttmp->GetBinContent(prochisttmp->GetMaximumBin())*1.1);
    secprochisttmp->SetMaximum(secprochisttmp->GetBinContent(secprochisttmp->GetMaximumBin())*1.1);


    //Normalizes 2D by X slice by calling NormalHist
    TH2D* proccorrhistnormtmp = NormalHist(proccorrhisttmp); //Doesn't really work i.e. entries not really spread out enough to warrant normalizing


    //Loops to label bins
    for(int jj=1;jj<=12;jj++)
    {
        prochisttmp->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
        proccorrhistnormtmp->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
    }

    for(int jj=1;jj<=8;jj++)
    {
        secprochisttmp->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        proccorrhistnormtmp->GetYaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    //More Styling
    prochisttmp->GetXaxis()->SetLabelSize(0.03);
    secprochisttmp->GetXaxis()->SetLabelSize(0.03);

    proccorrhistnormtmp->GetXaxis()->SetLabelSize(0.03);
    proccorrhistnormtmp->GetXaxis()->SetTitleOffset(0.8);
    //proccorrhistnormtmp->GetYaxis()->SetLabelOffset(0.8); //Unnecessary?

    //Generate name from cuts
    TString name = Cut2Name(allcut);


    //Print histos
    proccorrhistnormtmp->Draw("COLZ");
    c1->Print(Form("outplot/%s/proccorr%s.png",outdir.Data(),name.Data()));
    prochisttmp->Draw();
    c1->Print(Form("outplot/%s/proc%s.png",outdir.Data(),name.Data()));
    secprochisttmp->Draw();
    c1->Print(Form("outplot/%s/secproc%s.png",outdir.Data(),name.Data()));



    //Assign input histograms to clones of temporary (passing by reference)
    proc = (TH1D*)prochisttmp->Clone(Form("proc%s",name.Data()));
    secproc = (TH1D*)secprochisttmp->Clone(Form("secproc%s",name.Data()));
    proccorr = (TH2D*)proccorrhistnormtmp->Clone(Form("proccorr%s",name.Data()));

    delete proccorrhistnormtmp;
    delete prochisttmp;
    delete secprochisttmp;
    delete proccorrhisttmp;

    return procnorm;

}


//Converts given cut string into a name for legends,filename,histnames etc.
TString PlotHelp::Cut2Name(const TString allcut)
{
    TString name("");


    if(allcut.Contains("&"))
    {
        Int_t pos = -99;
        TString tmpcut = allcut;
        while(tmpcut.Contains("<"))
        {
            pos = tmpcut.First("<");
            name += "_"+tmpcut(pos+1,2);
            tmpcut.Remove(pos,1);
        }

    }
    else if(allcut.Contains("<"))
    {
            name += "_"+allcut(allcut.First("<")+1,2);
    }
    else
    {
        name += "Default";
    }

    return name; //name should contain just the numbers of the cut (tbd. remove . from decimals to avoid file problems)
}




//Takes array of histograms and plots them on same pad 
void PlotHelp::StackProc(TH1D* histo[],TCanvas* c1,TList* lout, const TString allcut[], 
                         const TString outdir,Double_t nor[],Int_t nhists,Bool_t kSec)
{
    TString name("");
    THStack* hs = 0x0;
    TLegend* lg = 0x0;

    //Split cases for secondary or primary process
    if(!kSec)
    {
        hs = new THStack("hsproc","Processes");
        lg = new TLegend(0.65,0.65,0.95,0.95);
    }
    else if(kSec)
    {
        hs = new THStack("hssecproc","Secondary Processes");
        lg = new TLegend(0.12,0.65,0.65,0.95);
    }
    
    lg->SetHeader("Cuts (Entries)");

    for(Int_t ii = 0; ii<nhists;ii++)
    {
        name = Cut2Name(allcut[ii]);
        hs->Add(histo[ii]);
        lg->AddEntry(histo[ii],Form("%s (%g)",name.Data(),nor[ii]/nor[0]));
    }

    style::ResetStyle(lg);

    c1->Clear();
    c1->SetGrid();
    hs->Draw("NOSTACK");
    lg->Draw();

    //Titles for X and Y axis must be assigned after Draw()
    hs->GetXaxis()->SetTitle("Process");
    hs->GetYaxis()->SetTitle("Frequency");
    if(kSec) hs->GetYaxis()->SetTitleOffset(1.2);
    c1->Modified();
    if(!kSec) c1->Print(Form("outplot/%s/procs.png",outdir.Data()));
    if(kSec)  c1->Print(Form("outplot/%s/secprocs.png",outdir.Data()));

    lout->Add(hs);

}



#endif

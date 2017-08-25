#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2D.h"
#include "TList.h"
#include "THStack.h"

#include "PlotHelp.h"

#include <cmath>

const char *capproclabel[8] = {"Absorption","Multi-Production", "Charge Exchange","Double Charge Exchange","Inelastic", 
    "Rangeout","Decay","Conversion"};

void ProcCmp()
{
    TList* lout = new TList;

    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");
    TH2D* rmaxhist = new TH2D("rmaxhist","Correlation between r_{max} and Process;Process;r_{max}/cm",8,0,8,19,5,100);
    TH2D* ovflowhist = new TH2D("ovflowhist","Correlation between Energy Overflow and Process;Process;Overflow/Mev",8,0,8,18,1,10.0);
    
    TCanvas* c1 = new TCanvas("c1","",900,600);
    tree->Draw("rmax:secondary_proc>>+rmaxhist");
    tree->Draw("enovflow:secondary_proc>>+ovflowhist");

    
    //Loops to label bins
    for(int jj=1;jj<=8;jj++)
    {
        rmaxhist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        ovflowhist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    rmaxhist->SetStats(0);
    ovflowhist->SetStats(0);

    TH2D* ovflownorm = PlotHelp::NormalHist(ovflowhist);
    TH2D* rmaxnorm = (TH2D*)rmaxhist->Clone();
    rmaxnorm->Scale(1/rmaxnorm->Integral());
    TH2D* rmaxcolnorm = PlotHelp::NormalHist(rmaxhist);

    lout->Add(ovflownorm);
    lout->Add(rmaxnorm);
    lout->Add(rmaxcolnorm);
    lout->Add(ovflowhist);

    TFile* fout = new TFile("outplot/ProCmp.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete c1;
}


void ERCorr()
{
    TList* lout = new TList;
    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");

    Double_t radius = -99.0;
    Double_t Esum = -99.0;

    Double_t Energybincut[300];
    Double_t rmax;
    Int_t    indmax;
    Int_t    secprocnum;

    tree->SetBranchAddress("Energybincut",&Energybincut);
    tree->SetBranchAddress("rmax",&rmax);
    tree->SetBranchAddress("indmax",&indmax);
    tree->SetBranchAddress("secondary_proc",&secprocnum);
   


    TH2D* erhists[8];
    TH2D* unerhists[8];
    TH2D* othhist = new TH2D("othhist","Others",20,0,100,25,0,150);

    //loop for othhist
    for(Long64_t jentry = 0;jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(secprocnum < 5)
        {
            Esum = 0;
            for(Int_t ii = 0 ;ii<indmax;ii++)
            {
                Esum += Energybincut[ii];
                radius = (ii+1)*5.0; //in cm
                othhist->Fill(rmax-radius,Esum);//Esum
            }

        }      
    }
    //loop over all processes
    for(Int_t iproc = 0; iproc<8; iproc++)
    {
        TH2D* tmph = new TH2D("tmph","",20,0,100,25,0,150);

        for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
        {
            tree->GetEntry(jentry);

            if(secprocnum == iproc)
            {
                Esum = 0;
                for(Int_t ii = 0 ;ii<indmax;ii++)
                {
                    Esum +=Energybincut[ii];
                    radius = (ii+1)*5.0; //in cm
                    tmph->Fill(rmax-radius,Esum);
                }

            }
        }

        tmph->SetXTitle("#Deltar = r_{max} - r  (cm)");
        tmph->SetYTitle("Energy Deposited Outside #Deltar  (MeV)");
        tmph->SetTitle(Form("Correlation between Energy Deposited and Relative Radius for %s ",capproclabel[iproc]));

       TH2D* normtmph = PlotHelp::NormalHist(tmph); 
       normtmph->SetStats(0);
       
       erhists[iproc] = (TH2D*)normtmph->Clone(secproclabel[iproc]);
       unerhists[iproc] = (TH2D*)tmph->Clone(Form("%s_un",secproclabel[iproc]));

       lout->Add(erhists[iproc]);
       lout->Add(unerhists[iproc]);

       delete tmph;
       delete normtmph;
    }

    lout->Add(othhist);

    TFile* fout = new TFile("outplot/Proc2DPlots.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
}


void HeadPlot()
{
    TList* lout = new TList;
    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");

    Double_t radius = -99.0;
    Double_t Esum = -99.0;

    Double_t Energybin[300];
    Int_t    indmax;
    Int_t    secprocnum;

    tree->SetBranchAddress("Energybin",&Energybin);
    tree->SetBranchAddress("indmax",&indmax);
    tree->SetBranchAddress("secondary_proc",&secprocnum);
   


    TH2D* erhists[8];
    TH2D* unerhists[8];
    TH2D* othhist = new TH2D("othhist","",20,0,100,25,0,150);


    //loop to fill othhist
    for(Long64_t jentry = 0;jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(secprocnum < 5)
        {
            Esum = 0;
            for(Int_t ii = 0;ii<100;ii++)
            {
                Esum +=Energybin[ii];
                radius = (ii)*5.0; //in cm
                othhist->Fill(radius,Esum);
            }

        }      
    }


    //loop over all processes
    for(Int_t iproc = 0; iproc<8; iproc++)
    {
        TH2D* tmph = new TH2D("tmph","",20,0,100,25,0,150);

        for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
        {
            tree->GetEntry(jentry);


            if(secprocnum == iproc)
            {
                Esum = 0;

                for(Int_t ii = 0;ii<100;ii++)
                {
                    Esum +=Energybin[ii];
                    radius = (ii)*5.0; //in cm 
                    tmph->Fill(radius,Esum);
                }
            }
        }

        tmph->SetXTitle("Radius from Track Endpoint (cm)");
        tmph->SetYTitle("Cumulative Energy Deposited (MeV)");
        tmph->SetTitle(Form("Correlation between Energy Deposited and Radius for %s ",capproclabel[iproc]));

       TH2D* normtmph = PlotHelp::NormalHist(tmph); 
       normtmph->SetStats(0);
       
       erhists[iproc] = (TH2D*)normtmph->Clone(secproclabel[iproc]);
       unerhists[iproc] = (TH2D*)tmph->Clone(Form("%s_un",secproclabel[iproc]));

       lout->Add(erhists[iproc]);
       lout->Add(unerhists[iproc]);

       delete tmph;
       delete normtmph;
    }

    lout->Add(othhist);

    TFile* fout = new TFile("outplot/Proc2DHeadPlots.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
}


void ProjPlot(const TString file, const TString outname)
{
    TList* lout = new TList;
    THStack* hs[20];
    TFile* fin = new TFile(Form("outplot/%s",file.Data()),"READ");
    TCanvas* c1 = new TCanvas("c1","",900,600);

    TH2D* rangehist =(TH2D*)fin->Get("rangeout_un");
    TH2D* othhist =(TH2D*)fin->Get("othhist");

//Loop over bins 
    for(Int_t jj = 0; jj < 20; jj++)
    {
        THStack* tmphs = 0x0;
        if(outname.Contains("tail")) tmphs =  new THStack(Form("radius_%g",jj*5.0),Form("Energy Deposited outside of  Relative Radius %g",jj*5.0));
        if(outname.Contains("head")) tmphs = new THStack(Form("radius_%g",jj*5.0),Form("Energy Deposited within Radius %g",(jj+1)*5.0));

        TLegend* leg = new TLegend(0.65,0.65,0.85,0.85);
        style::ResetStyle(leg);

        TH1D* tmprange = (TH1D*)rangehist->ProjectionY("range_py",jj+1,jj+1);
        TH1D* tmpoth = (TH1D*)othhist->ProjectionY("oth_py",jj+1,jj+1);

        tmprange->Scale(1/tmprange->Integral());
        tmpoth->Scale(1/tmpoth->Integral());

        tmprange->SetLineWidth(2);
        tmpoth->SetLineWidth(2);

        tmprange->SetLineColor(kRed);
        tmpoth->SetLineColor(kBlack);

        tmphs->Add(tmprange);
        tmphs->Add(tmpoth);

        leg->AddEntry(tmprange,"Rangeout");
        leg->AddEntry(tmpoth,"Others");

        tmphs->Draw("NOSTACK");
        leg->Draw();
        tmphs->GetXaxis()->SetTitle("Cumulative Energy Deposited (MeV)");
        tmphs->GetYaxis()->SetTitle("Frequency");
        c1->Modified();

        if(outname.Contains("tail"))
        {
            c1->Print(Form("outplot/%s/radius_%g.png",outname.Data(),jj*5.0));
            hs[jj] = (THStack*)tmphs->Clone(Form("radius_%g",jj*5));
        }


        if(outname.Contains("head"))
        {
            c1->Print(Form("outplot/%s/radius_%g.png",outname.Data(),(jj+1)*5.0));
            hs[jj] = (THStack*)tmphs->Clone(Form("radius_%g",(jj+1)*5));
        }


        lout->Add(hs[jj]);

        delete tmphs;
        delete leg;
        //delete tmpharr;
    }

    TFile* fout = new TFile(Form("outplot/%s/RadProc.root",outname.Data()),"RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
    delete c1;
}


void ZeroPlot()
{
    const TString trkrcut("(ecal==0)&&(ecal10==0)&&(ecal15==0)&&(ecal20==0)&&(ecal25==0)&&(ecal30==0)");

    TList* lout = new TList;
    TTree* tree = PlotHelp::InputFile("Pion_Output_MINERvA.root","tout");
    TH1D* tmphfail = new TH1D("zerohist","Effect of Zero Suppression",8,0,8);
    TH1D* tmphall = new TH1D("allhist","",8,0,8);
    TH1D* tmphpass = new TH1D("passhist","",8,0,8);

    TLegend* leg = new TLegend(0.70,0.65,0.90,0.90);

    TCanvas* c1 = new TCanvas("c1","",900,600);


    tree->Draw("secondary_proc>>+zerohist","(Energybin[0]==0)&&"+trkrcut);
    tree->Draw("secondary_proc>>+allhist",trkrcut);
    tree->Draw("secondary_proc>>+passhist","(Energybin[0]>0)&&"+trkrcut);


    Double_t noall = tmphall->GetEntries();
    Double_t nopass = tmphpass->GetEntries();
    Double_t nofail = tmphfail->GetEntries();

    tmphall->Scale(1/noall);
    tmphpass->Scale(1/nopass);
    tmphfail->Scale(1/nofail);

    tmphfail->SetStats(0);
    tmphpass->SetStats(0);
    tmphall->SetStats(0);

    tmphfail->SetXTitle("Process");
    tmphfail->SetYTitle("Events");
    tmphall->SetXTitle("Process");
    tmphall->SetYTitle("Events");
    tmphpass->SetXTitle("Process");
    tmphpass->SetYTitle("Events");

    tmphfail->SetLineWidth(2);
    tmphall->SetLineWidth(2);
    tmphpass->SetLineWidth(2);


    tmphall->SetLineColor(kBlack);
    tmphfail->SetLineColor(kRed);
    tmphpass->SetLineColor(kGreen);

    leg->AddEntry(tmphall,Form("All ( %.3g )",noall/noall));
    leg->AddEntry(tmphpass,Form("Pass ( %.3g )",nopass/noall));
    leg->AddEntry(tmphfail,Form("Fail ( %.3g )",nofail/noall));

    style::ResetStyle(leg);

    for(int jj=1;jj<=8;jj++)
    {
        tmphfail->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        tmphall->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        tmphpass->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    tmphfail->Draw();
    tmphpass->Draw("SAME");
    tmphall->Draw("SAME");
    leg->Draw();

    c1->Print("outplot/zerosup.png");

    lout->Add(tmphfail);
    lout->Add(tmphall);
    lout->Add(tmphpass);


    TFile* fout = new TFile("outplot/Zeroplot.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
    delete c1;
}


void CutPlot()
{
    TList* lout = new TList;

    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");

    TH1D* histo = new TH1D("histo","Processes before  and after cut (Area Normalized)",8,0,8);
    TH1D* hist = new TH1D("hist","Processes after r_{max} cut (Area Normalized)",8,0,8);

    Int_t secprocnum, indmax;
    Double_t Energybincut[150];


    tree->SetBranchAddress("indmax",&indmax);
    tree->SetBranchAddress("Energybincut",&Energybincut);
    tree->SetBranchAddress("secondary_proc",&secprocnum);

    for(Long64_t jentry = 0;jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        if(secprocnum != 6)  ////Filter out decay events
        {
        histo->Fill(secprocnum);

        Double_t Esum = 0.0;
        for(Int_t ii = 0;ii<indmax;ii++)
        {
            Esum +=Energybincut[ii];
        }

        if(Esum<35) hist->Fill(secprocnum);
        }
    }


    Double_t scal = hist->Integral();
    Double_t scalo = histo->Integral();


    hist->Scale(1/scal);
    histo->Scale(1/scalo);

    histo->SetLineColor(kBlue);
    hist->SetLineColor(kRed);

    hist->SetStats(0);
    histo->SetStats(0);

    TLegend* leg = new TLegend(0.65,0.65,0.85,0.90);
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

void AllProjPlot(const TString file, const TString outname)
{
    TList* lout = new TList;
    THStack* hs[20];
    TFile* fin = new TFile(Form("outplot/%s",file.Data()),"READ");
    TCanvas* c1 = new TCanvas("c1","",900,600);


//Loop over bins 
    for(Int_t jj = 0; jj < 20; jj++)
    {
        THStack* tmphs = 0x0;
        if(outname.Contains("tail")) tmphs =  new THStack(Form("radius_%g",jj*5.0),Form("Energy Deposited outside of Relative Radius %g",jj*5.0));
        if(outname.Contains("head")) tmphs = new THStack(Form("radius_%g",jj*5.0),Form("Energy Deposited within Radius %g",(jj+1)*5.0));

        TLegend* leg = new TLegend(0.65,0.65,0.85,0.85);
        style::ResetStyle(leg);

        TH1D* tmpharr[8];

        for(Int_t ii = 0; ii <7; ii++)
        {
            TH2D* tmp2h =(TH2D*)fin->Get(Form("%s_un",secproclabel[ii]));
            TH1D* tmph = (TH1D*)tmp2h->ProjectionY(Form("%s_py",secproclabel[ii]),jj+1,jj+1);

            if(tmph->Integral() > 0) tmph->Scale(1/tmph->Integral());

            tmph->SetLineWidth(2);

            tmph->SetLineColor(ColorArray[ii]);

            tmpharr[ii] = (TH1D*)tmph->Clone();

            tmphs->Add(tmpharr[ii]);

            leg->AddEntry(tmpharr[ii],secproclabel[ii]);

            delete tmph;
            delete tmp2h;
        }

        tmphs->Draw("NOSTACK");
        leg->Draw();
        tmphs->GetXaxis()->SetTitle("Cumulative Energy Deposited (MeV)");
        tmphs->GetYaxis()->SetTitle("Frequency");
        c1->Modified();

        if(outname.Contains("tail"))
        {
            c1->Print(Form("outplot/%s/radius_%g.png",outname.Data(),jj*5.0));
            hs[jj] = (THStack*)tmphs->Clone(Form("radius_%g",jj*5));
        }


        if(outname.Contains("head"))
        {
            c1->Print(Form("outplot/%s/radius_%g.png",outname.Data(),(jj+1)*5.0));
            hs[jj] = (THStack*)tmphs->Clone(Form("radius_%g",(jj+1)*5));
        }


        lout->Add(hs[jj]);

        delete tmphs;
        delete leg;
        //delete tmpharr;
    }

    TFile* fout = new TFile(Form("outplot/%s/RadProc.root",outname.Data()),"RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();

    delete lout;
    delete fout;
    delete c1;
}

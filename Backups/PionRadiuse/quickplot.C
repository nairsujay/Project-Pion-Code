#include "PlotHelp.h"
#include "style.h"

void quickplot()
{
    TList* lout = new TList;

    TTree* tree = PlotHelp::InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");

    TH1D* histo = new TH1D("histo","Processes before  and after cut (Area Normalized)",8,0,8);
    TH1D* hist = new TH1D("hist","Processes after r_{max} cut (Area Normalized)",8,0,8);

    TCanvas* c1 = new TCanvas("c1","",900,600);

    tree->Draw("secondary_proc>>+hist","rmax<10");
    tree->Draw("secondary_proc>>+histo");

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
    leg->AddEntry(hist,Form("r_{max} < 10 ( %g )",scal/scalo));


    for(int jj=1;jj<=8;jj++)
    {
        hist->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        histo->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    }

    style::ResetStyle(leg);

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

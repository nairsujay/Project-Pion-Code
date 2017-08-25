#include <TH1.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>



void doCmp2()
{
    TH1D *Ainelhist[5];
    TH1D *Arangehist[5];
    TH1D *Binelhist[5];
    TH1D *Brangehist[5];
    TH1D* diffhist[5];
    TH1D* ratiohist[5];
    TList* lout = new TList;
    
//Debug related
    //Int_t ine=0; //number of bins
    //Int_t ren=0; //number of bins
    //Double_t subint=0;
    //Double_t divint=0;
    //Get Histograms


    for(Int_t ii=0; ii<5;ii++)
    {
        Ainelhist[ii] = 0x0;
        Arangehist[ii] = 0x0;


        TH1D* ineltmp =(TH1D*) gDirectory->Get(Form("ineltrkr%d",(ii+2)*5));
        TH1D* rangetmp =(TH1D*) gDirectory->Get(Form("rangeouttrkr%d",(ii+2)*5));


        if(!ineltmp)
        {
            printf("no inel \n");
            exit(1);
        } 

        if(!rangetmp)
        {
            printf("no rangeout \n");
            exit(1);
        } 



        Ainelhist[ii] = (TH1D*)ineltmp->Clone(Form("Aineltrkr%d",(ii+2)*5));
        Arangehist[ii] = (TH1D*)rangetmp->Clone(Form("Arangetrkr%d",(ii+2)*5));
        Binelhist[ii] = (TH1D*)ineltmp->Clone(Form("Bineltrkr%d",(ii+2)*5));
        Brangehist[ii] = (TH1D*)rangetmp->Clone(Form("Brangetrkr%d",(ii+2)*5));


        //Normalizes A histograms by area
        Ainelhist[ii]->Scale(1/(Ainelhist[ii]->Integral()));
        Arangehist[ii]->Scale(1/(Arangehist[ii]->Integral()));
        

        //Normalizes B histograms by max bin
        Binelhist[ii]->Scale(1/(Binelhist[ii]->GetBinContent(Binelhist[ii]->GetMaximumBin())));
        Brangehist[ii]->Scale(1/(Brangehist[ii]->GetBinContent(Brangehist[ii]->GetMaximumBin())));

        lout->Add(Ainelhist[ii]);
        lout->Add(Arangehist[ii]);
        lout->Add(Binelhist[ii]);
        lout->Add(Brangehist[ii]);


        delete ineltmp;
        delete rangetmp;
    }



    //Manipulate Histograms
    for(Int_t jj=0; jj<5;jj++)
    {
        Int_t radius = (jj+2)*5;
        diffhist[jj] = 0x0;
        ratiohist[jj] = 0x0;


        TH1D* difftmp = new TH1D("inelsubtmp","",60,0,160);
        TH1D* ratiotmp = new TH1D("rangedivtmp","",60,0,160);
        TCanvas* c1 = new TCanvas("c1","",600,400);


       difftmp->Add(Binelhist[jj],Brangehist[jj],1,-1);
       ratiotmp->Divide(Ainelhist[jj],Arangehist[jj],1,1);


       difftmp->SetTitle(Form("ineltrkr%d-rangetrkr%d;Energy;Occurences",radius,radius));
       ratiotmp->SetTitle(Form("ineltrkr%d/rangetrkr%d;Energy;Occurences",radius,radius));


       difftmp->SetMaximum(difftmp->GetBinContent(difftmp->GetMaximumBin())*1.1);
       ratiotmp->SetMaximum(ratiotmp->GetBinContent(ratiotmp->GetMaximumBin())*1.1);
       

       diffhist[jj] = (TH1D*)difftmp->Clone(Form("difftrkr%d",radius));
       ratiohist[jj] = (TH1D*)ratiotmp->Clone(Form("ratiotrkr%d",radius));


       diffhist[jj]->Draw("HIST");
       c1->Print(Form("outplot/difftrkr_%d.png",radius));
       ratiohist[jj]->Draw("HIST");
       c1->Print(Form("outplot/ratiotrkr_%d.png",radius));

       lout->Add(diffhist[jj]);
       lout->Add(ratiohist[jj]);


       delete difftmp;
       delete ratiotmp;
       delete c1;


    }


    //File to save 
    TFile* fout = new TFile("outplot/CmpPlot.root","UPDATE");
    lout->Write();
    fout->Save();
    fout->Close();
    delete fout;


}




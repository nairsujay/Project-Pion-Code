#include <TH1.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>


void doCmp()
{
    TH1D *inelhist[5];
    TH1D *rangehist[5];
    TH1D* inelsub[5];
    TH1D* rangesub[5];
    TH1D* ineldiv[5];
    TH1D* rangediv[5];
    TList* lout = new TList;
    
//Debug related
    Int_t ine=0; //number of bins
    Int_t ren=0; //number of bins
    Double_t ransubint=0;
    Double_t randivint=0;
    //Get Histograms


    for(Int_t ii=0; ii<5;ii++)
    {
        inelhist[ii] = 0x0;
        rangehist[ii] = 0x0;


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

        inelhist[ii] = (TH1D*)ineltmp->Clone(Form("ineltrkr%d",ii));
        rangehist[ii] = (TH1D*)rangetmp->Clone(Form("rangetrkr%d",ii));

        
        lout->Add(inelhist[ii]);
        lout->Add(rangehist[ii]);


        delete ineltmp;
        delete rangetmp;
    }



    //Manipulate Histograms
    for(Int_t jj=0; jj<5;jj++)
    {
        inelsub[jj] = 0x0;
        rangesub[jj] = 0x0;
        ineldiv[jj] = 0x0;
        rangediv[jj] = 0x0;


        ine = inelhist[jj]->GetSize();
        ren = rangehist[jj]->GetSize();
        

        printf("%d\t%d\n",ine,ren);


        TH1D* inelsubtmp = new TH1D("inelsubtmp","",60,0,160);
        TH1D* rangesubtmp = new TH1D("rangesubtmp","",60,0,160);
        TH1D* ineldivtmp = new TH1D("ineldivtmp","",60,0,160);
        TH1D* rangedivtmp = new TH1D("rangedivtmp","",60,0,160);
        TCanvas* c1 = new TCanvas("c1","",600,400);


        if(jj==0)
        {
            inelsubtmp->Add(inelhist[4],inelhist[0],1.0,-1.0);
            rangesubtmp->Add(rangehist[4],rangehist[0],1.0,-1.0);
            ineldivtmp->Divide(inelhist[4],inelhist[0],1.0,1.0);
            rangedivtmp->Divide(rangehist[4],rangehist[0],1.0,1.0);


            inelsub[jj]=(TH1D*)inelsubtmp->Clone(Form("ineltrkrsub%d%d",30,10));
            rangesub[jj]=(TH1D*)rangesubtmp->Clone(Form("rangetrkrsub%d%d",30,10));
            ineldiv[jj]=(TH1D*)ineldivtmp->Clone(Form("ineltrkrdiv%d%d",30,10));
            rangediv[jj]=(TH1D*)rangedivtmp->Clone(Form("rangetrkrdiv%d%d",30,10));


            inelsub[jj]->SetTitle(Form("ineltrkr%d-ineltrkr%d;Energy;Events",30,10));
            rangesub[jj]->SetTitle(Form("rangetrkr%d-rangetrkr%d;Energy;Events",30,10));
            ineldiv[jj]->SetTitle(Form("ineltrkr%d/ineltrkr%d;Energy;Events",30,10));
            rangediv[jj]->SetTitle(Form("rangetrkr%d/rangetrkr%d;Energy;Events",30,10));


            inelsub[jj]->SetAxisRange(5.,160.,"X");
            inelsub[jj]->SetAxisRange(-2.,2.,"Y");
            rangesub[jj]->SetAxisRange(5.,160.,"X");
            rangesub[jj]->SetAxisRange(-2.,2.,"Y");
            ineldiv[jj]->SetAxisRange(5.,160.,"X");
            ineldiv[jj]->SetAxisRange(-2.,2.,"Y");
            rangediv[jj]->SetAxisRange(5.,160.,"X");
            rangediv[jj]->SetAxisRange(-2.,2.,"Y");


            ransubint=rangesub[jj]->Integral();
            randivint=rangediv[jj]->Integral();
            printf("%g\t%g\n",ransubint,randivint);

            
            inelsub[jj]->Scale(1.0/(inelsub[jj]->Integral()));
            rangesub[jj]->Scale(1.0/(rangesub[jj]->Integral()));
            ineldiv[jj]->Scale(1.0/(ineldiv[jj]->Integral()));
            rangediv[jj]->Scale(1.0/(rangediv[jj]->Integral()));


            inelsub[jj]->SetMaximum(inelsub[jj]->GetBinContent(inelsub[jj]->GetMaximumBin())*1.1);
            rangesub[jj]->SetMaximum(rangesub[jj]->GetBinContent(rangesub[jj]->GetMaximumBin())*1.1);
            ineldiv[jj]->SetMaximum(ineldiv[jj]->GetBinContent(ineldiv[jj]->GetMaximumBin())*1.1);
            rangediv[jj]->SetMaximum(rangediv[jj]->GetBinContent(rangediv[jj]->GetMaximumBin())*1.1);


            inelsub[jj]->Draw("HIST");
            c1->Print(Form("outplot/inelsub_%d_%d.png",30,10)); 
            rangesub[jj]->Draw("HIST");
            c1->Print(Form("outplot/rangesub_%d_%d.png",30,10)); 
            ineldiv[jj]->Draw("HIST");
            c1->Print(Form("outplot/ineldiv_%d_%d.png",30,10));
            rangediv[jj]->Draw("HIST");
            c1->Print(Form("outplot/rangediv_%d_%d.png",30,10));


            delete inelsubtmp;
            delete rangesubtmp;
            delete ineldivtmp;
            delete rangedivtmp;
            delete c1;
            lout->Add(inelsub[jj]);
            lout->Add(rangesub[jj]);
            lout->Add(ineldiv[jj]);
            lout->Add(rangediv[jj]);


            continue;

        }

        inelsubtmp->Add(inelhist[jj],inelhist[jj-1],1.0,-1.0);
        rangesubtmp->Add(rangehist[jj],rangehist[jj-1],1.0,-1.0);
        ineldivtmp->Divide(inelhist[jj],inelhist[jj-1],1.0,1.0);
        rangedivtmp->Divide(rangehist[jj],rangehist[jj-1],1.0,1.0);


        inelsub[jj]=(TH1D*)inelsubtmp->Clone(Form("ineltrkrsub%d%d",(jj+2)*5,(jj+1)*5));
        rangesub[jj]=(TH1D*)rangesubtmp->Clone(Form("rangetrkrsub%d%d",(jj+2)*5,(jj+1)*5));
        ineldiv[jj]=(TH1D*)ineldivtmp->Clone(Form("ineltrkrdiv%d%d",(jj+2)*5,(jj+1)*5));
        rangediv[jj]=(TH1D*)rangedivtmp->Clone(Form("rangetrkrdiv%d%d",(jj+2)*5,(jj+1)*5));


        //Sets histo Style and normalizes
        inelsub[jj]->SetTitle(Form("ineltrkr%d-ineltrkr%d;Energy;Events",(jj+2)*5,(jj+1)*5));
        rangesub[jj]->SetTitle(Form("rangetrkr%d-rangetrkr%d;Energy;Events",(jj+2)*5,(jj+1)*5));
        ineldiv[jj]->SetTitle(Form("ineltrkr%d/ineltrkr%d;Energy;Events",(jj+2)*5,(jj+1)*5));
        rangediv[jj]->SetTitle(Form("rangetrkr%d/rangetrkr%d;Energy;Events",(jj+2)*5,(jj+1)*5));

        inelsub[jj]->SetAxisRange(5.,160.,"X");
        inelsub[jj]->SetAxisRange(-2.,2.,"Y");
        rangesub[jj]->SetAxisRange(5.,160.,"X");
        rangesub[jj]->SetAxisRange(-2.,2.,"Y");
        ineldiv[jj]->SetAxisRange(5.,160.,"X");
        ineldiv[jj]->SetAxisRange(-2.,2.,"Y");
        rangediv[jj]->SetAxisRange(5.,160.,"X");
        rangediv[jj]->SetAxisRange(-2.,2.,"Y");

        ransubint=rangesub[jj]->Integral();
        randivint=rangediv[jj]->Integral();
        printf("%g\t%g\n",ransubint,rangediv);

        inelsub[jj]->Scale(1.0/(inelsub[jj]->Integral()));
        rangesub[jj]->Scale(1.0/(rangesub[jj]->Integral()));
        ineldiv[jj]->Scale(1.0/(ineldiv[jj]->Integral()));
        rangediv[jj]->Scale(1.0/(rangediv[jj]->Integral()));

        inelsub[jj]->SetMaximum(inelsub[jj]->GetBinContent(inelsub[jj]->GetMaximumBin())*1.1);
        rangesub[jj]->SetMaximum(rangesub[jj]->GetBinContent(rangesub[jj]->GetMaximumBin())*1.1);
        ineldiv[jj]->SetMaximum(ineldiv[jj]->GetBinContent(ineldiv[jj]->GetMaximumBin())*1.1);
        rangediv[jj]->SetMaximum(rangediv[jj]->GetBinContent(rangediv[jj]->GetMaximumBin())*1.1);



        inelsub[jj]->Draw("HIST");
        c1->Print(Form("outplot/inelsub_%d_%d.png",(jj+2)*5,(jj+1)*5)); 
        rangesub[jj]->Draw("HIST");
        c1->Print(Form("outplot/rangesub_%d_%d.png",(jj+2)*5,(jj+1)*5)); 
        ineldiv[jj]->Draw("HIST");
        c1->Print(Form("outplot/ineldiv_%d_%d.png",(jj+2)*5,(jj+1)*5));
        rangediv[jj]->Draw("HIST");
        c1->Print(Form("outplot/rangediv_%d_%d.png",(jj+2)*5,(jj+1)*5));

        lout->Add(inelsub[jj]);
        lout->Add(rangesub[jj]);
        lout->Add(ineldiv[jj]);
        lout->Add(rangediv[jj]);


        delete inelsubtmp;
        delete rangesubtmp;
        delete ineldivtmp;
        delete rangedivtmp;
        delete c1;
       

    }


    //File to save 
    TFile* fout = new TFile("outplot/CmpPlot.root","RECREATE");
    lout->Write();
    fout->Save();
    fout->Close();
    delete fout;


}




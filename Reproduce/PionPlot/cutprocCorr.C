
void cutprocCorr()
{

    Int_t proc_num = -99;
    Int_t secproc_num = -99;
    //Double_t trkr30 =-999.0;
    //Double_t threshold = 0.0;
    const TString allcut("(trkr30 < 50)&&(trkr25 < 40)&&(trkr20 < 30)&&(trkr15 < 20)&&(trkr10 < 10)");

   // TFile::Open("Pion_Output_MINERvA.root");

    TTree* tree = (TTree*)gDirectory->Get("tout");
    //tree->SetBranchAddress("proc",&proc_num);
    //tree->SetBranchAddress("secondary_proc",&secproc_num);
    //tree->SetBranchAddress("trkr30",&trkr30);


    
    
    const char *proclabel[12] = {"hadElastic","NeutronInelastic","ProtonInelastic",
                                 "PionMinusInelastic","PionPlusInelastic",
                                 "CHIPS","Decay","conv", "compt", 
                                 "Primary","Other","Rangeout"};
    const char *secproclabel[8] = {"abs","multiprod", "cex","dbl_cex","inel", 
                                    "rangeout","decay","conv"};


        TH1I* prochist = new TH1I("prochist","Process Frequencies",12,0,12);
        TH1I* secprochist = new TH1I("secprochist","Secondary-Process Frequencies",8,0,8);
        TH2I* proccorrhist = new TH2I("proccorrhist","Process and Secondary Process",12,0,12,8,0,8);

    //main loop for filling histogram
    //for(int ii=0;ii<tree->GetEntries();ii++)
    //{
        //tree->GetEntry(ii);
        ////cout<<proc_num<<"\t"<<secproc_num<<endl;


        //if((trkr30 < 50)&&(trkr25 < 40)&&(trkr20 < 30)&&(trkr15 < 20)&&(trkr10 < 10))
        //{
            //prochist->Fill(proc_num);
            //secprochist->Fill(secproc_num);
            //proccorrhist->Fill(proc_num,secproc_num);
        //}
    //}

    TCanvas* c2 = new TCanvas("c2","");
    tree->Draw("proc>>+prochist",allcut);
    tree->Draw("secondary_proc>>+secprochist",allcut);
    tree->Draw("secondary_proc:proc>>+proccorrhist",allcut);



    //Styling
    //prochist->SetStats(0);
    //secprochist->SetStats(0);
    proccorrhist->SetStats(0);

    prochist->SetFillColor(4);
    prochist->SetXTitle("Process");
    prochist->SetYTitle("Events");

    secprochist->SetFillColor(4);
    secprochist->SetXTitle("Secondary Process");
    secprochist->SetYTitle("Events");

    prochist->SetXTitle("Process");
    prochist->SetYTitle("Secondary Process");



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

    TFile* fout = new TFile("outplot/cutproccorr/cutproccorr.root","RECREATE");
    prochist->Write();
    secprochist->Write();
    proccorrhist->Write();
    fout->Close();
    TCanvas* c1= new TCanvas("c1","proc");
    c1->SetGrid();
    prochist->Draw();
    c1->Print("outplot/cutproccorr/cutproc.png");
    secprochist->Draw();
    c1->Print("outplot/cutproccorr/cutsecproc.png");
    proccorrhist->Draw("COLZ");
    c1->Print("outplot/cutproccorr/cutproccorr.png");
}



//Double_t ProcPlot(TTree* tree,TCanvas* c1, const TString allcut,const TString outdir,
              //TH1D* &proc, TH1D* &secproc, TH1D* &proccorr) 
//{
    //const char *proclabel[12] = {"hadElastic","NeutronInelastic","ProtonInelastic",
        //"PionMinusInelastic","PionPlusInelastic",
        //"CHIPS","Decay","conv", "compt", 
        //"Primary","Other","Rangeout"};
    //const char *secproclabel[8] = {"abs","multiprod", "cex","dbl_cex","inel", 
        //"rangeout","decay","conv"};

    //const Color_t ColourArray[4] = {kBlack,kBlue,kRed,kGreen};


    //Double_t procnorm = -99;
    //Double_t checknorm = -99;


    //TH1D* prochisttmp = new TH1D("prochisttmp","Process Frequencies",12,0,12);
    //TH1D* secprochisttmp = new TH1D("secprochisttmp","Secondary-Process Frequencies",8,0,8);
    //TH2D* proccorrhisttmp = new TH2D("proccorrhisttmp",
                                     //"Process and Secondary Process Correlation",
                                     //12,0,12,8,0,8);

    
    //tree->Draw("proc>>+prochisttmp",allcut);
    //tree->Draw("secondary_proc>>+secprochisttmp",allcut);
    //tree->Draw("secondary_proc:proc>>+proccorrhisttmp",allcut);


    ////Styling

    //style::ResetStyle(prochisttmp);
    //style::ResetStyle(secprochisttmp);
    //style::ResetStyle(proccorrhisttmp);

    //prochisttmp->SetStats(0);
    //secprochisttmp->SetStats(0);
    //proccorrhisttmp->SetStats(0);

    //prochisttmp->GetXaxis()->SetTickLength(0);
    //secprochisttmp->GetXaxis()->SetTickLength(0);

    //prochisttmp->SetLineWidth(2);
    //secprochisttmp->SetLineWidth(2);

    ////prochisttmp->SetXTitle("Process");
    ////prochisttmp->SetYTitle("Events");

    ////secprochisttmp->SetXTitle("Secondary");
    ////secprochisttmp->SetYTitle("Events");
    
    //procnorm = prochisttmp->GetEntries();
    //checknorm = secprochisttmp->GetEntries();


////Exit if secproc and proc have different entries
    //if(procnorm != checknorm){printf("\n\ndifferent norm\n\n");exit(1);};



    ////normalize 1D histos by area and style them further
    //prochisttmp->Scale(1.0/procnorm[ii]);
    //secprochisttmp->Scale(1.0/secnorm[ii]);

    //prochisttmp->SetMaximum(prochisttmp->GetBinContent(prochisttmp->GetMaximumBin())*1.1);
    //secprochisttmp->SetMaximum(secprochisttmp->GetBinContent(secprochisttmp->GetMaximumBin())*1.1);


    ////Normalizes 2D by X slice
    //TH2D* proccorrhistnormtmp = NormalHist(proccorrhisttmp); //Doesn't really work...


    ////Loops to label bins
    //for(int jj=1;jj<=12;jj++)
    //{
        //prochisttmp->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
        //proccorrhistnormtmp->GetXaxis()->SetBinLabel(jj,proclabel[jj-1]);
    //}

    //for(int jj=1;jj<=8;jj++)
    //{
        //secprochisttmp->GetXaxis()->SetBinLabel(jj,secproclabel[jj-1]);
        //proccorrhistnormtmp->GetYaxis()->SetBinLabel(jj,secproclabel[jj-1]);
    //}

    //prochisttmp->GetXaxis()->SetLabelSize(0.03);
    //secprochisttmp->GetXaxis()->SetLabelSize(0.03);

    //proccorrhistnormtmp->GetXaxis()->SetLabelSize(0.03);
    //proccorrhistnormtmp->GetXaxis()->SetTitleOffset(0.8);
    ////proccorrhistnormtmp->GetYaxis()->SetLabelOffset(0.8); //Unnecessary?

    //prochisttmp->SetLineColor(ColourArray[ii]);
    //secprochisttmp->SetLineColor(ColourArray[ii]);

    ////Print tmp histos
    //proccorrhistnormtmp->Draw("COLZ");
    //c1->Print(Form("outplot/%s/proccorr_%s.png",outdir,allcut));
    //prochisttmp->Draw();
    //c1->Print(Form("outplot/%s/proc_%s.png",outdir,allcut));
    //secprochisttmp->Draw();
    //c1->Print(Form("outplot/%s/secproc_%s.png",outdir,allcut));


    
    //proc = (TH1D*)prochisttmp->Clone(Form("proc%s",allcut));
    //secproc = (TH1D*)secprochisttmp->Clone(Form("secproc%s",allcut));
    //proccorr = (TH2D*)proccorrhistnormtmp->Clone(Form("proccorr%s",allcut));


    //return procnorm;

//}

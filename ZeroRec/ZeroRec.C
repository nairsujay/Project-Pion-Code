#include "TTree.h"
#include "TFile.h"

#include <vector>



TTree* InputFile(const TString fname,const TString tname)
{
    TFile* fin = new TFile(fname,"READ");
    TTree* tree =(TTree*)fin->Get(tname);

    return tree;
}


void ZeroRec()
{
    TFile* file = new TFile("Pion_Output_Recur.root","UPDATE");
    TTree* tree = (TTree*)file->Get("tout");

    Double_t Energybin[200];

    tree->SetBranchAddress("Energybin",Energybin);

    Int_t nrec = 0; //number of recursions required

    Int_t modsize = 0; //size of vector

    std::vector<double> Energymodbin; //Will eventually contain energybins after being recursed through

    TBranch* nrecbranch = tree->Branch("nrec",&nrec,"nrec/I");
    TBranch* energymodbranch = tree->Branch("Energymodbin",&Energymodbin); 
    TBranch* modsizebranch = tree->Branch("modsize",&modsize,"modsize/I");

   //Main for loop loops through tree
   for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
   {
       tree->GetEntry(jentry);

       Energymodbin.clear();
       nrec = 0;
       modsize = 0;

       //skips unnecessary assignments etc
       if(Energybin[0] > 0)
       {
           Energymodbin.assign(Energybin,Energybin+200); 
           nrec = 0;
           modsize = Energymodbin.size();        
           nrecbranch->Fill();
           modsizebranch->Fill();
           energymodbranch->Fill();
           continue;
       } 

       Double_t esum = 0;
       for(Int_t ii = 0; ii<200;ii++)
       {
           esum+=Energybin[ii];
       }
       

       //Special case where string of bins is all zeros --- Prevents segfaults during pop_back
       if(esum == 0)
       {
           Energymodbin.assign(Energybin,Energybin+200); 
           nrec = -99;
           modsize = 0;        
           nrecbranch->Fill();
           modsizebranch->Fill();
           energymodbranch->Fill();
           continue;
       }
         

       std::vector<double> Energyrevbin; //Energybin with contents reversed for pop_back
   
    
       //Loops over Energybin array to populate reverse vector
       for(Int_t ii = 199;ii>=0;ii+=-1)
       {
           Energyrevbin.push_back(Energybin[ii]);
       }



       //Loops over Energyrevbin removing empty bins at end
       while(Energyrevbin.back() == 0)
       {
           Energyrevbin.pop_back();
           nrec++;
       }


       //Restores original order
       for(Int_t ii = Energyrevbin.size()-1;ii>=0;ii+=-1)
       {
           Energymodbin.push_back(Energyrevbin[ii]);
       }

       modsize = Energymodbin.size();

       nrecbranch->Fill();
       modsizebranch->Fill();
       energymodbranch->Fill();

       //Energyrevbin should now contain the energies but recursed until first bin is non-zero

   }

   tree->Write("",TObject::kOverwrite);
   file->Save();
   file->Close();
   
   delete file;
    
}

void ZeroFilter()
{
    TTree* intree = InputFile("Pion_Output_Recur.root","tout");

    //Filtered tree
    TFile* fout = new TFile("Filtered_Pion_Output_trkr_nozero.root","RECREATE");
    TTree* filteree =(TTree*)intree->CopyTree("(nrec>=0)&&(truth_Trkr)"); 
    
    //Selects entries for which at least one energy bin is filled and event is in trkr
    
    filteree->Write();
    fout->Close();
    delete fout;
}



void PrintCheck()
{
    //TTree* tree = InputFile("Filtered_Pion_Output_trkr_nozero.root","tout");
    TTree* tree = InputFile("Pion_Output_Recur.root","tout");

    std::vector<double> *Energyarr = new std::vector<double>;
    Int_t nrec;

    tree->SetBranchAddress("Energymodbin",&Energyarr);
    tree->SetBranchAddress("nrec",&nrec);

    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);

        if(nrec <0)
        {

            //if((*Energyarr)[0] == 0){fprintf(stderr,"\n\n \t FAILED \t \n\n");exit(1);};
            printf("\n");
            printf(" %d ",Energyarr->size());
            //for(Int_t ii =0; ii<modsize;ii++)
            //{
                //printf(" %g ",(*Energyarr)[ii]);
            //}
        }
    }

    fprintf(stderr,"\n\n \t Success \t \n\n");
}




Int_t FindRadMax(Double_t binarray[], Int_t arsize)
{
    Int_t radmax = -99;

    for(Int_t ii = 0; ii<(arsize-1); ii++)
    {
        if(binarray[ii] == 0 && binarray[ii+1] == 0){radmax = ii; break;};
    }

    if(radmax == -99){fprintf(stderr,"\n\n infinite \n\n");exit(1);};

    return radmax;

}



void RadAna()
{
    TFile* file = new TFile("Filtered_Pion_Output_trkr_nozero.root","UPDATE");
    TTree* tree = (TTree*)file->Get("tout");


    Int_t modsize;
    std::vector<double> *Energyarr = new std::vector<double>;

    tree->SetBranchAddress("modsize",&modsize);
    tree->SetBranchAddress("Energymodbin",&Energyarr);

    Double_t rmax, enovflow;
    Int_t indmax;

    Double_t Energybincut[200] = {0.0};

    TBranch *rmaxbranch = tree->Branch("rmax",&rmax,"rmax/D");
    TBranch *indmaxbranch = tree->Branch("indmax",&indmax,"indmax/I"); 
    TBranch *enovflowbranch = tree->Branch("enovflow",&enovflow,"enovflow/D");
    TBranch *energycutbranch = tree->Branch("Energybincut",&Energybincut,"Energybincut[200]/D");



    for(Long64_t jentry = 0; jentry<tree->GetEntries();jentry++)
    {
        tree->GetEntry(jentry);
        //for(Int_t ii = 0;ii<radiusrange;ii++){printf(" %d ",Radiusbin[ii]);};
        indmax = FindRadMax(Energyarr->data(),modsize);
        rmax = indmax*5.0; //in cm 
        //printf("%g \n",rmax);

        enovflow = 0;

        for(Int_t ii = 0; ii<200;ii++){Energybincut[ii] = 0;}; //Initialize energybincut array

        for(Int_t ii = 0; ii<indmax; ii++){Energybincut[ii] = (*Energyarr)[ii];};
        for(Int_t ii = indmax;ii<modsize;ii++){enovflow+=(*Energyarr)[ii];};/*if(Energybin[ii] > 0) printf("\n abcdefg %llu %d %g %g  \n",jentry,ii,enovflow,Energybin[ii]);};*/


        //test 
        //if(enovflow > 0) exit(1);


        //Prints arrays to log 
        printf("\n");
        for(Int_t ii = 0; ii<indmax; ii++){printf(" %g ", Energybincut[ii]);};
        printf(" \t %g \t\t",enovflow);
        for(Int_t ii = 0; ii<indmax+6; ii++){printf(" %g ",(*Energyarr)[ii]);};
        printf("\t\t");



        rmaxbranch->Fill();
        indmaxbranch->Fill();
        energycutbranch->Fill();
        enovflowbranch->Fill();
}


tree->Write("",TObject::kOverwrite);
file->Save();
file->Close();

delete file;

}


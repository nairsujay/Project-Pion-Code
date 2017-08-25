#include "TTree.h"
#include "TFile.h"

#include <vector>
#include <cmath>


Double_t*  NeighFindRad( std::vector<Double_t> *z,std::vector<Double_t> *t, Int_t arrsize)
{
    Double_t* neighrad = new Double_t[arrsize];
    for(Int_t jj = 0; jj<arrsize;jj++)
    {
        //Position of cluster
        Double_t z0 = (*z)[jj];
        Double_t t0 = (*t)[jj];
        Double_t rad = 1e6;

        for(Int_t ii = 0; ii<arrsize;ii++)
        {
            if(ii != jj)
            {
                Double_t tmprad = sqrt(pow(((*z)[ii]-z0),2.0)+pow(((*t)[ii]-t0),2.0));
                if (tmprad<rad) rad = tmprad;
            }
        }

        neighrad[jj] = rad/10.0; //convert to cm
    }
        
    return neighrad;
}


Int_t* NeighFindInd(std::vector<Double_t> *z,std::vector<Double_t> *t, Int_t arrsize)
{
    Int_t* neighind = new Int_t[arrsize];

    for(Int_t jj = 0; jj<arrsize;jj++)
    {
        //Position of cluster
        Double_t z0 = (*z)[jj];
        Double_t t0 = (*t)[jj];
        Double_t rad = 1e6;
        Int_t ind = -99;

        for(Int_t ii = 0; ii<arrsize;ii++)
        {
            if(ii != jj)
            {
                Double_t tmprad = sqrt(pow(((*z)[ii]-z0),2.0)+pow(((*t)[ii]-t0),2.0));
                if (tmprad<rad){rad = tmprad;ind = ii;};
            }
        }

        neighind[jj] = ind;
    }

    return neighind;
}

void NeighWrite()
{
    TFile* file = new TFile("Pion_Output_Recur.root","UPDATE");
    TTree* tree = (TTree*)file->Get("tout");

    std::vector<Double_t> *Zarray = new std::vector<Double_t>;
    std::vector<Double_t> *Tarray = new std::vector<Double_t>;

    Int_t arrsize;
    std::vector<Int_t> neighind;
    std::vector<Double_t> neighrad;

    tree->SetBranchAddress("Zarray",&Zarray);
    tree->SetBranchAddress("Tarray",&Tarray);
    tree->SetBranchAddress("arrsize",&arrsize);

    TBranch* neighradbranch = tree->Branch("neighrad",&neighrad);
    TBranch* neighindbranch = tree->Branch("neighind",&neighind);

    for(Long64_t jentry = 0; jentry<tree->GetEntries(); jentry++)
    {
        tree->GetEntry(jentry);

        Int_t* tmpneighind = NeighFindInd(Zarray,Tarray,arrsize);
        Double_t* tmpneighrad = NeighFindRad(Zarray,Tarray,arrsize);

        neighind.assign(tmpneighind,tmpneighind+arrsize);
        neighrad.assign(tmpneighrad,tmpneighrad+arrsize);

        neighradbranch->Fill();
        neighindbranch->Fill();

        printf("\n");
        for(Int_t ii = 0;ii<arrsize; ii++){printf(" %g ",neighrad[ii]);};
    }

    tree->Write("",TObject::kOverwrite);
    file->Save();
    file->Close();

    delete file;
}

#define dEdxAnaTool_cxx
#include "dEdxAnaTool.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
#include <fstream>

void dEdxAnaTool::SetMap()
{
__mass_map[2212] = 938.2720;
__mass_map[-211] = 139.5702;
__mass_map[+211] = 139.5702;


__proc_map[0]  = "hadElastic";
__proc_map[1]  = "NeutronInelastic";
__proc_map[2]  = "ProtonInelastic";
__proc_map[3]  = "PionMinusInelastic";
__proc_map[4]  = "PionPlusInelastic";
__proc_map[5]  = "CHIPSNuclearCaptureAtRest";
__proc_map[6]  = "Decay";
__proc_map[7]  = "conv";
__proc_map[8]  = "compt";
__proc_map[10] = "Primary";
__proc_map[20] = "Other";

__procnum_map["hadElastic"] = 0;
__procnum_map["NeutronInelastic"]  = 1;
__procnum_map["ProtonInelastic"]  = 2;
__procnum_map["PionMinusInelastic"]  = 3;
__procnum_map["PionPlusInelastic"]  = 4;
__procnum_map["CHIPSNuclearCaptureAtRest"]  = 5;
__procnum_map["Decay"]  = 6;
__procnum_map["conv"]  = 7;
__procnum_map["compt"]  = 8;
__procnum_map["Primary"] = 9;
__procnum_map["Other"] = 10;
__procnum_map["Rangeout"] = 11;

__secprocnum_map["abs"] = 0;
__secprocnum_map["multiprod"]  = 1;
__secprocnum_map["cex"]  = 2;
__secprocnum_map["dbl_cex"]  = 3;
__secprocnum_map["inel"]  = 4;
__secprocnum_map["rangeout"]  = 5;
__secprocnum_map["decay"]  = 6;
__secprocnum_map["conv"]  = 7;

}

TChain * dEdxAnaTool::InputFiles(const TString file, const TString tr)
{
  TChain *ch=new TChain(tr);

  if(file.Contains(".root"))
    ch->Add(file);
  else{
    ifstream fin(file);
    if(!fin){
      printf("NeutrinoTools::InputFiles file not found \n%s\n\n",file.Data()); exit(1);
    }

    TString buff;
    while(fin.good()){
      fin>>buff;
      if(buff!=""){
        ch->Add(buff);
      }
    }
  }

  //const Int_t ent = ch->GetEntries(); //takes infinity time!!                                                                                                                                                            
  printf("\t%d trees!\n",ch->GetNtrees());

  return ch;
}



void dEdxAnaTool::Loop(const TString filelist)
{
//   In a ROOT session, you can do:
//      Root > .L dEdxAnaTool.C
//      Root > dEdxAnaTool t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

    
    
    Init(InputFiles(filelist,"dEdxAnaTool"));

//Check    
    if (fChain == 0) return;

//Initializes the map;
    SetMap();



//Define Branches for tree and other variables

//-------------------------------------------------------
//-------------------------------------------------------
    Int_t    pdg         = 0;
    Int_t    proc_num    = -99;
    Int_t    secproc_num    = -99;
    Double_t evis_trkr   = 0.0;
    Double_t evis_ecal   = 0.0;
    Double_t evis_trkr10 = 0.0;
    Double_t evis_ecal10 = 0.0;
    Double_t evis_trkr15 = 0.0;
    Double_t evis_ecal15 = 0.0;
    Double_t evis_trkr20 = 0.0;
    Double_t evis_ecal20 = 0.0;
    Double_t evis_trkr25 = 0.0;
    Double_t evis_ecal25 = 0.0;
    Double_t evis_trkr30 = 0.0;
    Double_t evis_ecal30 = 0.0;
    Double_t ecalo       = 0.0;
    Double_t E_E         = 0.0; 
    Double_t kinetic     = 0.0;
    Double_t Eint        = 0.0;
    std::string proc("Unknown");
    std::string secondary_proc("unknown");

//-------------------------------------------------------
//-------------------------------------------------------

//Define Tree and branches

//-------------------------------------------------------
//-------------------------------------------------------
    TFile* outfile=new TFile("Pion_Output_MINERvA.root","recreate");
    TTree* tout =new TTree("tout","Output Tree");
    tout->Branch("pdg",&pdg,"PDG/I");
    tout->Branch("proc",&proc_num,"proc/I");
    tout->Branch("secondary_proc",&secproc_num,"secondary_proc/I");
    tout->Branch("trkr",&evis_trkr,"Energy/D");
    tout->Branch("trkr10",&evis_trkr10,"Energy/D");
    tout->Branch("trkr15",&evis_trkr15,"Energy/D");
    tout->Branch("trkr20",&evis_trkr20,"Energy/D");
    tout->Branch("trkr25",&evis_trkr25,"Energy/D");
    tout->Branch("trkr30",&evis_trkr30,"Energy/D");
    tout->Branch("ecal",&evis_ecal,"Energy/D");
    tout->Branch("ecal10",&evis_ecal10,"Energy/D");
    tout->Branch("ecal15",&evis_ecal15,"Energy/D");
    tout->Branch("ecal20",&evis_ecal20,"Energy/D");
    tout->Branch("ecal25",&evis_ecal25,"Energy/D");
    tout->Branch("ecal30",&evis_ecal30,"Energy/D");
    tout->Branch("ecalo",&ecalo,"Energy/D");
    tout->Branch("E",&E_E,"Energy/D");
    tout->Branch("kinetic",&kinetic,"KineticEnergy/D");
    tout->Branch("Eint",&Eint,"KineticEnergy/D");
//-------------------------------------------------------
 
//-------------------------------------------------------

//Start of Loop through Entries
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      evis_trkr   = 0.0;
      evis_ecal   = 0.0;
      evis_trkr10 = 0.0;
      evis_ecal10 = 0.0;
      evis_trkr15 = 0.0;
      evis_ecal15 = 0.0;
      evis_trkr20 = 0.0;
      evis_ecal20 = 0.0;
      evis_trkr25 = 0.0;
      evis_ecal25 = 0.0;
      evis_trkr30 = 0.0;
      evis_ecal30 = 0.0;
      ecalo       = 0.0;
      E_E         = 0.0;
      kinetic     = 0.0;
      Eint        = 0.0;



      if (recovtx_ntrack > 0) {
          // remove pions with large z residual, I guess mostly elastic scattered
          // pions with the second segment not tracked
          if (abs(prong_stop_z-detmc_traj_zf[0]) > 100.0) continue;
      }

      pdg = mc_FSPartPDG[0];
      E_E = mc_FSPartE[0]; //renamed E to E_E to stop warning from ACLiC
      kinetic = E_E - __mass_map[pdg];
      Eint = detmc_traj_preEf[0] - __mass_map[pdg];


      int npip = 0;
      int npi0 = 0;
      int npim = 0;
      std::vector<int> secondaries;
      for (int i = 0; i < detmc_ntrajectory2; ++i) {
          if (detmc_traj_mother[i] != 1) continue;
          if (detmc_traj_proc[i] == 0) continue; // hadElastic;

          secondaries.push_back(detmc_traj_pdg[i]);

          if (detmc_traj_pdg[i] == 211) ++npip;
          if (detmc_traj_pdg[i] ==-211) ++npim;
          if (detmc_traj_pdg[i] == 111) ++npi0;

//#define print_it 
//#ifdef print_it
          //printf("\t %3d %3d %15d %30s %10.2f %10.2f %10.2f %10.2f\n",
                  //detmc_traj_id[i],
                  //detmc_traj_mother[i],
                  //detmc_traj_pdg[i],
                  //__proc_map[detmc_traj_proc[i]].c_str(),
                  //detmc_traj_E0[i],
                  //detmc_traj_x0[i],
                  //detmc_traj_y0[i],
                  //detmc_traj_z0[i]);
//#endif

          proc = __proc_map[detmc_traj_proc[i]];
      }

      if      (npim + npi0 + npip == 0) secondary_proc = "abs";
      else if (npim + npi0 + npip > 1) secondary_proc = "multiprod";
      else if (npim == 0 && npi0 == 1 && npip == 0) secondary_proc = "cex";
      else if (pdg == 211 && (npim == 1 && npi0 == 0 && npip == 0)) secondary_proc = "dbl_cex";
      else if (pdg ==-211 && (npim == 0 && npi0 == 0 && npip == 1)) secondary_proc = "dbl_cex";
      else if (pdg == 211 && (npim == 0 && npi0 == 0 && npip == 1)) secondary_proc = "inel";
      else if (pdg ==-211 && (npim == 1 && npi0 == 0 && npip == 0)) secondary_proc = "inel";

      if (pdg == 2212 && proc == "Unknown") proc = "Rangeout";

      if (pdg ==  211 && Eint < 1.0) proc = "Rangeout";
      if (pdg == -211 && proc == "CHIPSNuclearCaptureAtRest") proc = "Rangeout";

      // Repeat the main process if the sub-process is not well-defined
      if (proc == "Rangeout") secondary_proc = "rangeout";
      if (proc == "Decay") secondary_proc = "decay";


      for (int i = 0; i < cluster_energy_vec_sz; ++i) {

          if (cluster_used_vec[i] == 1) continue;

          double energy = cluster_energy_vec[i];
          double z = cluster_z_vec[i];
          double t = cluster_x_vec[i];

          double z0 = prong_stop_z;
          double t0 = prong_stop_x;
          if (cluster_view_vec[i] == 2) t0 = prong_stop_u;
          if (cluster_view_vec[i] == 3) t0 = prong_stop_v;

          double radius = sqrt(pow(z-z0,2) + pow(t-t0,2))/10.0; // in cm

          if   (z < 8590.0) evis_trkr += energy;
          else evis_ecal += energy;

          if (z < 8590.0) {
              if (radius < 10.0) evis_trkr10 += energy;
              if (radius < 15.0) evis_trkr15 += energy;
              if (radius < 20.0) evis_trkr20 += energy;
              if (radius < 25.0) evis_trkr25 += energy;
              if (radius < 30.0) evis_trkr30 += energy;

          } 
          else {
              if (radius < 10.0) evis_ecal10 += energy;
              if (radius < 15.0) evis_ecal15 += energy;
              if (radius < 20.0) evis_ecal20 += energy;
              if (radius < 25.0) evis_ecal25 += energy;
              if (radius < 30.0) evis_ecal30 += energy;
          }

      }

      ecalo = evis_trkr20/0.8 + evis_ecal20/0.45;
      if (pdg == 211 && ecalo > 30.0 && secondary_proc == "rangeout") {
          for (int i = 0; i < cluster_energy_vec_sz; ++i) {

              double energy = cluster_energy_vec[i];
              double z = cluster_z_vec[i];
              double t = cluster_x_vec[i];

              //printf("\t %4d %4d %2d %10.3f %10.3f %10.3f\n", i,
                      //cluster_used_vec[i], cluster_view_vec[i],
                      //z, t, energy);
          }

      }
      if(pdg == 211 && evis_ecal == 0 && evis_trkr10 != 0)
      {
        printf("%10llu %d %10.2f %10.2f %30s %10s %4d %10.2f %5d %5d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
               jentry, pdg, kinetic, Eint, proc.c_str(), secondary_proc.c_str(), recovtx_ntrack,
               llr_score_proton, llr_offset_proton, llr_offset_pion,
               prong_hypo_pion_E-__mass_map[pdg],
               evis_trkr,
               evis_ecal,
               evis_trkr10,
               evis_ecal10,
               evis_trkr15,
               evis_ecal15,
               evis_trkr20,
               evis_ecal20,
               evis_trkr25,
               evis_ecal25,
               evis_trkr30,
               evis_ecal30,
               track_stop_zvec[0]);
      }

      proc_num = __procnum_map[proc];
      secproc_num = __secprocnum_map[secondary_proc];

      tout->Fill();
   }
   outfile->cd();
   tout->Write();
   outfile->Save();
   outfile->Close();
}

#define selector_dedx_calo_cxx
// The class definition in selector_dedx_calo.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("selector_dedx_calo.C")
// Root > T->Process("selector_dedx_calo.C","some options")
// Root > T->Process("selector_dedx_calo.C+")
//
#include <iostream>
#include <cmath>


#include "selector_dedx_calo.h"
#include <TH2.h>
#include <TStyle.h>
#include <TRandom1.h>
#include "Event.h"
#include <TTree.h>
#include <TFile.h>

int readEntry = 0;
TRandom1 generator(123456789);


void selector_dedx_calo::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   
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


}

void selector_dedx_calo::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void selector_dedx_calo::Loop()
{
    TFile outfile("Test.root","RECREATE");
    TTree* tout =new TTree("tout","Output Tree");
    Event* levent=new Event(); 

    const Int_t ent = fChain->GetEntries();
    for(int ii=0; ii<ent; ii++){
        cout<<"Entry "<<ii<<endl;
        Process(ii);
    }
    outfile.Close();

}


Bool_t selector_dedx_calo::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either selector_dedx_calo::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  //fChain->GetTree("dEdxAnaTool")->GetEntry(entry);
  fChain->GetEntry(entry);

    //int ran_val = generator.Uniform(0,100);

    //if (ran_val > 10) return false;
    
    ++readEntry;

    if (recovtx_ntrack > 0) {
            // remove pions with large z residual, I guess mostly elastic scattered
            // pions with the second segment not tracked
        if (abs(prong_stop_z-detmc_traj_zf[0]) > 100.0) return false;
    }
    
    int pdg = mc_FSPartPDG[0];
    double E = mc_FSPartE[0];
    double kinetic = E - __mass_map[pdg];
    double Eint = detmc_traj_preEf[0] - __mass_map[pdg];

    std::string proc("Unknown");
    

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

#define print_it 
#ifdef print_it
        printf("\t %3d %3d %15d %30s %10.2f %10.2f %10.2f %10.2f\n",
               detmc_traj_id[i],
               detmc_traj_mother[i],
               detmc_traj_pdg[i],
               __proc_map[detmc_traj_proc[i]].c_str(),
               detmc_traj_E0[i],
               detmc_traj_x0[i],
               detmc_traj_y0[i],
               detmc_traj_z0[i]);
#endif
        proc = __proc_map[detmc_traj_proc[i]];
    }

  std::string secondary_proc("unknown");

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

    double evis_trkr = 0.0;
    double evis_ecal = 0.0;

    double evis_trkr10 = 0.0;
    double evis_ecal10 = 0.0;
    double evis_trkr15 = 0.0;
    double evis_ecal15 = 0.0;
    double evis_trkr20 = 0.0;
    double evis_ecal20 = 0.0;
    double evis_trkr25 = 0.0;
    double evis_ecal25 = 0.0;
    double evis_trkr30 = 0.0;
    double evis_ecal30 = 0.0;

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

        } else {
            if (radius < 10.0) evis_ecal10 += energy;
            if (radius < 15.0) evis_ecal15 += energy;
            if (radius < 20.0) evis_ecal20 += energy;
            if (radius < 25.0) evis_ecal25 += energy;
            if (radius < 30.0) evis_ecal30 += energy;
        }
        
    }

    double ecalo = evis_trkr20/0.8 + evis_ecal20/0.45;
    if (pdg == 211 && ecalo > 30.0 && secondary_proc == "rangeout") {
        for (int i = 0; i < cluster_energy_vec_sz; ++i) {

            double energy = cluster_energy_vec[i];
            double z = cluster_z_vec[i];
            double t = cluster_x_vec[i];

            printf("\t %4d %4d %2d %10.3f %10.3f %10.3f\n", i,
                   cluster_used_vec[i], cluster_view_vec[i],
                   z, t, energy);
        }
        
    }
    
    printf("%10d %d %10.2f %10.2f %30s %10s %4d %10.2f %5d %5d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
           readEntry, pdg, kinetic, Eint, proc.c_str(), secondary_proc.c_str(), recovtx_ntrack,
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


    if (pdg ==-211 && (npim == 1 && npi0 == 0 && npip == 0))
    {
        Event* event=new Event();
        event->SetTrkr(evis_trkr,evis_trkr10,evis_trkr15,evis_trkr20,evis_trkr25,evis_trkr30);
        
        event->SetEcal(evis_ecal,evis_ecal10,evis_ecal15,evis_ecal20,evis_ecal25,evis_ecal30);
        tout->Fill();
        delete event;
    }


   // exit(1); //What's the purpose of this???
   return kTRUE;


}

void selector_dedx_calo::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void selector_dedx_calo::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Apr 16 23:30:30 2017 by ROOT version 5.34/32
// from TChain dEdxAnaTool/
//////////////////////////////////////////////////////////

#ifndef selector_dedx_calo_h
#define selector_dedx_calo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <map>

using std::vector;


class selector_dedx_calo : public TSelector {
public :
    std::map<int, double> __mass_map;
    std::map<int, std::string> __proc_map;

TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Double_t        eventID;
   Int_t           physEvtNum;
   Int_t           n_hyps;
   Int_t           processType;
   Int_t           primaryPart;
   Int_t           n_slices;
   Int_t           slice_numbers[1];   //[n_slices]
   Int_t           shared_slice;
   Double_t        vtx[4];
   Double_t        vtxErr[4];
   Double_t        E[4];
   Bool_t          found_truth;
   Bool_t          phys_front_activity;
   Bool_t          phys_energy_in_road_upstream_is_rockmuon_consistent;
   Bool_t          rock_muons_removed;
   Bool_t          minos_track_match;
   Bool_t          minos_stub_match;
   Bool_t          unknown_helicity;
   Bool_t          minos_track_inside_partial_plane;
   Bool_t          prim_vtx_has_misassigned_track_direction;
   Bool_t          prim_vtx_has_broken_track;
   Int_t           broken_track_most_us_plane;
   Int_t           detmc_ntrajectory;
   Int_t           detmc_ntrajectory2;
   Int_t           llr_offset_pion;
   Int_t           llr_offset_proton;
   Int_t           minos_trk_end_plane;
   Int_t           minos_trk_is_contained;
   Int_t           minos_trk_is_ok;
   Int_t           minos_trk_quality;
   Int_t           minos_trk_used_curvature;
   Int_t           minos_trk_used_range;
   Int_t           ncluster_all;
   Int_t           ncluster_used;
   Int_t           nprong;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_vertex_is_fiducial;
   Int_t           prong_ntrack;
   Int_t           recovtx_ntrack;
   Int_t           recovtx_ntrack2;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        llr_score_pion;
   Double_t        llr_score_proton;
   Double_t        minos_trk_eqp;
   Double_t        minos_trk_fit_pass;
   Double_t        minos_trk_p;
   Double_t        minos_trk_p_curvature;
   Double_t        minos_trk_p_range;
   Double_t        minos_trk_qp;
   Double_t        muon_E;
   Double_t        muon_dp;
   Double_t        muon_p;
   Double_t        muon_phi;
   Double_t        muon_px;
   Double_t        muon_py;
   Double_t        muon_pz;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        muon_thetax;
   Double_t        muon_thetay;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        pim_gl;
   Double_t        pim_kinetic;
   Double_t        pim_phi;
   Double_t        pim_pienergy;
   Double_t        pim_pimom;
   Double_t        pim_pl;
   Double_t        pim_pl_area;
   Double_t        pim_theta;
   Double_t        pim_thetax;
   Double_t        pim_thetay;
   Double_t        pip_gl;
   Double_t        pip_kinetic;
   Double_t        pip_phi;
   Double_t        pip_pienergy;
   Double_t        pip_pimom;
   Double_t        pip_pl;
   Double_t        pip_pl_area;
   Double_t        pip_theta;
   Double_t        pip_thetax;
   Double_t        pip_thetay;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        prong_hypo_pion_E;
   Double_t        prong_hypo_pion_charge;
   Double_t        prong_hypo_pion_mass;
   Double_t        prong_hypo_pion_p;
   Double_t        prong_hypo_pion_px;
   Double_t        prong_hypo_pion_py;
   Double_t        prong_hypo_pion_pz;
   Double_t        prong_hypo_pion_score;
   Double_t        prong_hypo_proton_E;
   Double_t        prong_hypo_proton_charge;
   Double_t        prong_hypo_proton_mass;
   Double_t        prong_hypo_proton_p;
   Double_t        prong_hypo_proton_px;
   Double_t        prong_hypo_proton_py;
   Double_t        prong_hypo_proton_pz;
   Double_t        prong_hypo_proton_score;
   Double_t        prong_length;
   Double_t        prong_length_area;
   Double_t        prong_start_u;
   Double_t        prong_start_v;
   Double_t        prong_start_x;
   Double_t        prong_start_y;
   Double_t        prong_start_z;
   Double_t        prong_stop_u;
   Double_t        prong_stop_v;
   Double_t        prong_stop_x;
   Double_t        prong_stop_y;
   Double_t        prong_stop_z;
   Int_t           cluster_used_vec_sz;
   Int_t           cluster_used_vec[113];   //[cluster_used_vec_sz]
   Int_t           cluster_view_vec_sz;
   Int_t           cluster_view_vec[113];   //[cluster_view_vec_sz]
   Int_t           detmc_traj_id_sz;
   Int_t           detmc_traj_id[40];   //[detmc_traj_id_sz]
   Int_t           detmc_traj_mother_sz;
   Int_t           detmc_traj_mother[40];   //[detmc_traj_mother_sz]
   Int_t           detmc_traj_pdg_sz;
   Int_t           detmc_traj_pdg[40];   //[detmc_traj_pdg_sz]
   Int_t           detmc_traj_proc_sz;
   Int_t           detmc_traj_proc[40];   //[detmc_traj_proc_sz]
   Int_t           detmc_traj_status_sz;
   Int_t           detmc_traj_status[40];   //[detmc_traj_status_sz]
   Int_t           track_module_vec_sz;
   Int_t           track_module_vec[110];   //[track_module_vec_sz]
   Int_t           track_ndigit_vec_sz;
   Int_t           track_ndigit_vec[110];   //[track_ndigit_vec_sz]
   Int_t           track_node_vec_sz;
   Int_t           track_node_vec[3];   //[track_node_vec_sz]
   Int_t           track_plane_vec_sz;
   Int_t           track_plane_vec[110];   //[track_plane_vec_sz]
   Int_t           track_truth_count_vec_sz;
   Int_t           track_truth_count_vec[3];   //[track_truth_count_vec_sz]
   Int_t           track_truth_pdg_vec1_sz;
   Int_t           track_truth_pdg_vec1[3];   //[track_truth_pdg_vec1_sz]
   Int_t           track_truth_pdg_vec2_sz;
   Int_t           track_truth_pdg_vec2[3];   //[track_truth_pdg_vec2_sz]
   Int_t           track_truth_pdg_vec3_sz;
   Int_t           track_truth_pdg_vec3[3];   //[track_truth_pdg_vec3_sz]
   Int_t           cluster_energy_vec_sz;
   Double_t        cluster_energy_vec[113];   //[cluster_energy_vec_sz]
   Int_t           cluster_time_vec_sz;
   Double_t        cluster_time_vec[113];   //[cluster_time_vec_sz]
   Int_t           cluster_x_vec_sz;
   Double_t        cluster_x_vec[113];   //[cluster_x_vec_sz]
   Int_t           cluster_z_vec_sz;
   Double_t        cluster_z_vec[113];   //[cluster_z_vec_sz]
   Int_t           detmc_traj_E0_sz;
   Double_t        detmc_traj_E0[40];   //[detmc_traj_E0_sz]
   Int_t           detmc_traj_Ef_sz;
   Double_t        detmc_traj_Ef[40];   //[detmc_traj_Ef_sz]
   Int_t           detmc_traj_preEf_sz;
   Double_t        detmc_traj_preEf[40];   //[detmc_traj_preEf_sz]
   Int_t           detmc_traj_prepxf_sz;
   Double_t        detmc_traj_prepxf[40];   //[detmc_traj_prepxf_sz]
   Int_t           detmc_traj_prepyf_sz;
   Double_t        detmc_traj_prepyf[40];   //[detmc_traj_prepyf_sz]
   Int_t           detmc_traj_prepzf_sz;
   Double_t        detmc_traj_prepzf[40];   //[detmc_traj_prepzf_sz]
   Int_t           detmc_traj_px0_sz;
   Double_t        detmc_traj_px0[40];   //[detmc_traj_px0_sz]
   Int_t           detmc_traj_pxf_sz;
   Double_t        detmc_traj_pxf[40];   //[detmc_traj_pxf_sz]
   Int_t           detmc_traj_py0_sz;
   Double_t        detmc_traj_py0[40];   //[detmc_traj_py0_sz]
   Int_t           detmc_traj_pyf_sz;
   Double_t        detmc_traj_pyf[40];   //[detmc_traj_pyf_sz]
   Int_t           detmc_traj_pz0_sz;
   Double_t        detmc_traj_pz0[40];   //[detmc_traj_pz0_sz]
   Int_t           detmc_traj_pzf_sz;
   Double_t        detmc_traj_pzf[40];   //[detmc_traj_pzf_sz]
   Int_t           detmc_traj_t0_sz;
   Double_t        detmc_traj_t0[40];   //[detmc_traj_t0_sz]
   Int_t           detmc_traj_tf_sz;
   Double_t        detmc_traj_tf[40];   //[detmc_traj_tf_sz]
   Int_t           detmc_traj_x0_sz;
   Double_t        detmc_traj_x0[40];   //[detmc_traj_x0_sz]
   Int_t           detmc_traj_xf_sz;
   Double_t        detmc_traj_xf[40];   //[detmc_traj_xf_sz]
   Int_t           detmc_traj_y0_sz;
   Double_t        detmc_traj_y0[40];   //[detmc_traj_y0_sz]
   Int_t           detmc_traj_yf_sz;
   Double_t        detmc_traj_yf[40];   //[detmc_traj_yf_sz]
   Int_t           detmc_traj_z0_sz;
   Double_t        detmc_traj_z0[40];   //[detmc_traj_z0_sz]
   Int_t           detmc_traj_zf_sz;
   Double_t        detmc_traj_zf[40];   //[detmc_traj_zf_sz]
   Int_t           recovtx_uvec_sz;
   Double_t        recovtx_uvec[6];   //[recovtx_uvec_sz]
   Int_t           recovtx_vvec_sz;
   Double_t        recovtx_vvec[6];   //[recovtx_vvec_sz]
   Int_t           recovtx_xvec_sz;
   Double_t        recovtx_xvec[6];   //[recovtx_xvec_sz]
   Int_t           recovtx_yvec_sz;
   Double_t        recovtx_yvec[6];   //[recovtx_yvec_sz]
   Int_t           recovtx_zvec_sz;
   Double_t        recovtx_zvec[6];   //[recovtx_zvec_sz]
   Int_t           track_cluster_energy_vec_sz;
   Double_t        track_cluster_energy_vec[11];   //[track_cluster_energy_vec_sz]
   Int_t           track_costheta_vec_sz;
   Double_t        track_costheta_vec[3];   //[track_costheta_vec_sz]
   Int_t           track_length_vec_sz;
   Double_t        track_length_vec[3];   //[track_length_vec_sz]
   Int_t           track_nodex_vec_sz;
   Double_t        track_nodex_vec[110];   //[track_nodex_vec_sz]
   Int_t           track_nodey_vec_sz;
   Double_t        track_nodey_vec[110];   //[track_nodey_vec_sz]
   Int_t           track_nodez_vec_sz;
   Double_t        track_nodez_vec[110];   //[track_nodez_vec_sz]
   Int_t           track_rev_cluster_energy_vec_sz;
   Double_t        track_rev_cluster_energy_vec[11];   //[track_rev_cluster_energy_vec_sz]
   Int_t           track_separation_vec_sz;
   Double_t        track_separation_vec[3];   //[track_separation_vec_sz]
   Int_t           track_slope_kx_vec_sz;
   Double_t        track_slope_kx_vec[110];   //[track_slope_kx_vec_sz]
   Int_t           track_slope_ky_vec_sz;
   Double_t        track_slope_ky_vec[110];   //[track_slope_ky_vec_sz]
   Int_t           track_start_uvec_sz;
   Double_t        track_start_uvec[3];   //[track_start_uvec_sz]
   Int_t           track_start_vvec_sz;
   Double_t        track_start_vvec[3];   //[track_start_vvec_sz]
   Int_t           track_start_xvec_sz;
   Double_t        track_start_xvec[3];   //[track_start_xvec_sz]
   Int_t           track_start_yvec_sz;
   Double_t        track_start_yvec[3];   //[track_start_yvec_sz]
   Int_t           track_start_zvec_sz;
   Double_t        track_start_zvec[3];   //[track_start_zvec_sz]
   Int_t           track_stop_uvec_sz;
   Double_t        track_stop_uvec[3];   //[track_stop_uvec_sz]
   Int_t           track_stop_vvec_sz;
   Double_t        track_stop_vvec[3];   //[track_stop_vvec_sz]
   Int_t           track_stop_xvec_sz;
   Double_t        track_stop_xvec[3];   //[track_stop_xvec_sz]
   Int_t           track_stop_yvec_sz;
   Double_t        track_stop_yvec[3];   //[track_stop_yvec_sz]
   Int_t           track_stop_zvec_sz;
   Double_t        track_stop_zvec[3];   //[track_stop_zvec_sz]
   Int_t           track_theta_vec_sz;
   Double_t        track_theta_vec[110];   //[track_theta_vec_sz]
   Int_t           track_thetax_vec_sz;
   Double_t        track_thetax_vec[110];   //[track_thetax_vec_sz]
   Int_t           track_thetay_vec_sz;
   Double_t        track_thetay_vec[110];   //[track_thetay_vec_sz]
   Int_t           track_truth_fraction_vec1_sz;
   Double_t        track_truth_fraction_vec1[3];   //[track_truth_fraction_vec1_sz]
   Int_t           track_truth_fraction_vec2_sz;
   Double_t        track_truth_fraction_vec2[3];   //[track_truth_fraction_vec2_sz]
   Int_t           track_truth_fraction_vec3_sz;
   Double_t        track_truth_fraction_vec3[3];   //[track_truth_fraction_vec3_sz]
   Int_t           track_truth_shared_vec_sz;
   Double_t        track_truth_shared_vec[3];   //[track_truth_shared_vec_sz]
   Bool_t          truth_has_physics_event;
   Bool_t          truth_is_fiducial;
   Bool_t          truth_pass_plausible;
   Bool_t          truth_is_signal;
   Int_t           truth_nneutron;
   Int_t           truth_nothermeson;
   Int_t           truth_npi0;
   Int_t           truth_npim;
   Int_t           truth_npip;
   Int_t           truth_nproton;
   Int_t           truth_one_pi0;
   Int_t           truth_one_pim;
   Int_t           truth_one_pip;
   Double_t        truth_fslepton_E;
   Double_t        truth_fslepton_P;
   Double_t        truth_fslepton_T;
   Double_t        truth_fslepton_phi;
   Double_t        truth_fslepton_theta;
   Double_t        truth_fslepton_theta_x;
   Double_t        truth_fslepton_theta_y;
   Int_t           dEdxAnaTool_nuFlavor;
   Int_t           dEdxAnaTool_nuHelicity;
   Int_t           dEdxAnaTool_intCurrent;
   Int_t           dEdxAnaTool_intType;
   Double_t        dEdxAnaTool_E;
   Double_t        dEdxAnaTool_Q2;
   Double_t        dEdxAnaTool_x;
   Double_t        dEdxAnaTool_y;
   Double_t        dEdxAnaTool_W;
   Double_t        dEdxAnaTool_score;
   Double_t        dEdxAnaTool_leptonE[4];
   Double_t        dEdxAnaTool_vtx[4];
   Int_t           ev_run;
   Int_t           ev_subrun;
   Int_t           ev_detector;
   Int_t           ev_triggerType;
   Int_t           ev_gate;
   Int_t           ev_global_gate;
   Int_t           ev_gps_time_sec;
   Int_t           ev_gps_time_usec;
   Int_t           mc_run;
   Int_t           mc_subrun;
   Int_t           mc_nInteractions;
   Int_t           mc_MIState;
   Double_t        mc_pot;
   Int_t           mc_beamConfig;
   Int_t           mc_processType;
   Int_t           mc_nthEvtInSpill;
   Int_t           mc_nthEvtInFile;
   Int_t           mc_intType;
   Int_t           mc_current;
   Int_t           mc_charm;
   Double_t        mc_weight;
   Double_t        mc_XSec;
   Double_t        mc_diffXSec;
   Int_t           mc_incoming;
   Double_t        mc_fluxDriverProb;
   Int_t           mc_targetNucleus;
   Int_t           mc_targetZ;
   Int_t           mc_targetA;
   Int_t           mc_targetNucleon;
   Int_t           mc_struckQuark;
   Int_t           mc_seaQuark;
   Int_t           mc_resID;
   Int_t           mc_primaryLepton;
   Double_t        mc_incomingE;
   Double_t        mc_Bjorkenx;
   Double_t        mc_Bjorkeny;
   Double_t        mc_Q2;
   Double_t        mc_nuT;
   Double_t        mc_w;
   Double_t        mc_vtx[4];
   Double_t        mc_incomingPartVec[4];
   Double_t        mc_initNucVec[4];
   Double_t        mc_primFSLepton[4];
   Int_t           mc_nFSPart;
   Double_t        mc_FSPartPx[1];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[1];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[1];   //[mc_nFSPart]
   Double_t        mc_FSPartE[1];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[1];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[1];   //[mc_er_nPart]
   Int_t           mc_er_status[1];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[1];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[1];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[1];   //[mc_er_nPart]
   Double_t        mc_er_Px[1];   //[mc_er_nPart]
   Double_t        mc_er_Py[1];   //[mc_er_nPart]
   Double_t        mc_er_Pz[1];   //[mc_er_nPart]
   Double_t        mc_er_E[1];   //[mc_er_nPart]
   Int_t           mc_er_FD[1];   //[mc_er_nPart]
   Int_t           mc_er_LD[1];   //[mc_er_nPart]
   Int_t           mc_er_mother[1];   //[mc_er_nPart]
   Int_t           mc_fr_nNuAncestorIDs;
   Int_t           mc_fr_nuAncestorIDs[1];   //[mc_fr_nNuAncestorIDs]
   Int_t           mc_fr_nuParentID;
   Int_t           mc_fr_decMode;
   Double_t        mc_fr_primProtonVtx[3];
   Double_t        mc_fr_primProtonP[4];
   Double_t        mc_fr_nuParentDecVtx[3];
   Double_t        mc_fr_nuParentProdVtx[3];
   Double_t        mc_fr_nuParentProdP[4];
   Double_t        mc_cvweight_total;
   Double_t        wgt;
   Double_t        mc_cvweight_totalFlux;
   Double_t        mc_cvweight_totalXsec;
   Double_t        mc_ppfx1_cvweight;
   Double_t        mc_hornCurrent_cvweight;
   Double_t        mc_gen1_cvweight_total;
   Double_t        gen1_wgt;
   Double_t        mc_gen1_cvweight_totalFlux;
   Double_t        mc_gen1_cvweight_NA49;
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[1];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_gen1_Flux_Tertiary_sz;
   Double_t        mc_wgt_gen1_Flux_Tertiary[1];   //[mc_wgt_gen1_Flux_Tertiary_sz]
   Int_t           mc_wgt_gen1_Flux_NA49_sz;
   Double_t        mc_wgt_gen1_Flux_NA49[1];   //[mc_wgt_gen1_Flux_NA49_sz]
   Int_t           mc_wgt_Norm_sz;
   Double_t        mc_wgt_Norm[1];   //[mc_wgt_Norm_sz]
   Int_t           mc_wgt_ppfx1_Total_sz;
   Double_t        mc_wgt_ppfx1_Total[1];   //[mc_wgt_ppfx1_Total_sz]
   Int_t           n_prongs;
   Int_t           prong_nParticles[3];   //[n_prongs]
   Double_t        prong_part_score[3];   //[n_prongs]
   Double_t        prong_part_mass[3];   //[n_prongs]
   Int_t           prong_part_charge[3];   //[n_prongs]
   Int_t           prong_part_pid[3];   //[n_prongs]
   vector<vector<double> > *prong_part_E;
   vector<vector<double> > *prong_part_pos;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_physEvtNum;   //!
   TBranch        *b_n_hyps;   //!
   TBranch        *b_processType;   //!
   TBranch        *b_primaryPart;   //!
   TBranch        *b_n_slices;   //!
   TBranch        *b_slice_numbers;   //!
   TBranch        *b_shared_slice;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxErr;   //!
   TBranch        *b_E;   //!
   TBranch        *b_found_truth;   //!
   TBranch        *b_phys_front_activity;   //!
   TBranch        *b_phys_energy_in_road_upstream_is_rockmuon_consistent;   //!
   TBranch        *b_rock_muons_removed;   //!
   TBranch        *b_minos_track_match;   //!
   TBranch        *b_minos_stub_match;   //!
   TBranch        *b_unknown_helicity;   //!
   TBranch        *b_minos_track_inside_partial_plane;   //!
   TBranch        *b_prim_vtx_has_misassigned_track_direction;   //!
   TBranch        *b_prim_vtx_has_broken_track;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_detmc_ntrajectory;   //!
   TBranch        *b_detmc_ntrajectory2;   //!
   TBranch        *b_llr_offset_pion;   //!
   TBranch        *b_llr_offset_proton;   //!
   TBranch        *b_minos_trk_end_plane;   //!
   TBranch        *b_minos_trk_is_contained;   //!
   TBranch        *b_minos_trk_is_ok;   //!
   TBranch        *b_minos_trk_quality;   //!
   TBranch        *b_minos_trk_used_curvature;   //!
   TBranch        *b_minos_trk_used_range;   //!
   TBranch        *b_ncluster_all;   //!
   TBranch        *b_ncluster_used;   //!
   TBranch        *b_nprong;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_prong_ntrack;   //!
   TBranch        *b_recovtx_ntrack;   //!
   TBranch        *b_recovtx_ntrack2;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_llr_score_pion;   //!
   TBranch        *b_llr_score_proton;   //!
   TBranch        *b_minos_trk_eqp;   //!
   TBranch        *b_minos_trk_fit_pass;   //!
   TBranch        *b_minos_trk_p;   //!
   TBranch        *b_minos_trk_p_curvature;   //!
   TBranch        *b_minos_trk_p_range;   //!
   TBranch        *b_minos_trk_qp;   //!
   TBranch        *b_muon_E;   //!
   TBranch        *b_muon_dp;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_muon_thetax;   //!
   TBranch        *b_muon_thetay;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_pim_gl;   //!
   TBranch        *b_pim_kinetic;   //!
   TBranch        *b_pim_phi;   //!
   TBranch        *b_pim_pienergy;   //!
   TBranch        *b_pim_pimom;   //!
   TBranch        *b_pim_pl;   //!
   TBranch        *b_pim_pl_area;   //!
   TBranch        *b_pim_theta;   //!
   TBranch        *b_pim_thetax;   //!
   TBranch        *b_pim_thetay;   //!
   TBranch        *b_pip_gl;   //!
   TBranch        *b_pip_kinetic;   //!
   TBranch        *b_pip_phi;   //!
   TBranch        *b_pip_pienergy;   //!
   TBranch        *b_pip_pimom;   //!
   TBranch        *b_pip_pl;   //!
   TBranch        *b_pip_pl_area;   //!
   TBranch        *b_pip_theta;   //!
   TBranch        *b_pip_thetax;   //!
   TBranch        *b_pip_thetay;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_prong_hypo_pion_E;   //!
   TBranch        *b_prong_hypo_pion_charge;   //!
   TBranch        *b_prong_hypo_pion_mass;   //!
   TBranch        *b_prong_hypo_pion_p;   //!
   TBranch        *b_prong_hypo_pion_px;   //!
   TBranch        *b_prong_hypo_pion_py;   //!
   TBranch        *b_prong_hypo_pion_pz;   //!
   TBranch        *b_prong_hypo_pion_score;   //!
   TBranch        *b_prong_hypo_proton_E;   //!
   TBranch        *b_prong_hypo_proton_charge;   //!
   TBranch        *b_prong_hypo_proton_mass;   //!
   TBranch        *b_prong_hypo_proton_p;   //!
   TBranch        *b_prong_hypo_proton_px;   //!
   TBranch        *b_prong_hypo_proton_py;   //!
   TBranch        *b_prong_hypo_proton_pz;   //!
   TBranch        *b_prong_hypo_proton_score;   //!
   TBranch        *b_prong_length;   //!
   TBranch        *b_prong_length_area;   //!
   TBranch        *b_prong_start_u;   //!
   TBranch        *b_prong_start_v;   //!
   TBranch        *b_prong_start_x;   //!
   TBranch        *b_prong_start_y;   //!
   TBranch        *b_prong_start_z;   //!
   TBranch        *b_prong_stop_u;   //!
   TBranch        *b_prong_stop_v;   //!
   TBranch        *b_prong_stop_x;   //!
   TBranch        *b_prong_stop_y;   //!
   TBranch        *b_prong_stop_z;   //!
   TBranch        *b_cluster_used_vec_sz;   //!
   TBranch        *b_cluster_used_vec;   //!
   TBranch        *b_cluster_view_vec_sz;   //!
   TBranch        *b_cluster_view_vec;   //!
   TBranch        *b_detmc_traj_id_sz;   //!
   TBranch        *b_detmc_traj_id;   //!
   TBranch        *b_detmc_traj_mother_sz;   //!
   TBranch        *b_detmc_traj_mother;   //!
   TBranch        *b_detmc_traj_pdg_sz;   //!
   TBranch        *b_detmc_traj_pdg;   //!
   TBranch        *b_detmc_traj_proc_sz;   //!
   TBranch        *b_detmc_traj_proc;   //!
   TBranch        *b_detmc_traj_status_sz;   //!
   TBranch        *b_detmc_traj_status;   //!
   TBranch        *b_track_module_vec_sz;   //!
   TBranch        *b_track_module_vec;   //!
   TBranch        *b_track_ndigit_vec_sz;   //!
   TBranch        *b_track_ndigit_vec;   //!
   TBranch        *b_track_node_vec_sz;   //!
   TBranch        *b_track_node_vec;   //!
   TBranch        *b_track_plane_vec_sz;   //!
   TBranch        *b_track_plane_vec;   //!
   TBranch        *b_track_truth_count_vec_sz;   //!
   TBranch        *b_track_truth_count_vec;   //!
   TBranch        *b_track_truth_pdg_vec1_sz;   //!
   TBranch        *b_track_truth_pdg_vec1;   //!
   TBranch        *b_track_truth_pdg_vec2_sz;   //!
   TBranch        *b_track_truth_pdg_vec2;   //!
   TBranch        *b_track_truth_pdg_vec3_sz;   //!
   TBranch        *b_track_truth_pdg_vec3;   //!
   TBranch        *b_cluster_energy_vec_sz;   //!
   TBranch        *b_cluster_energy_vec;   //!
   TBranch        *b_cluster_time_vec_sz;   //!
   TBranch        *b_cluster_time_vec;   //!
   TBranch        *b_cluster_x_vec_sz;   //!
   TBranch        *b_cluster_x_vec;   //!
   TBranch        *b_cluster_z_vec_sz;   //!
   TBranch        *b_cluster_z_vec;   //!
   TBranch        *b_detmc_traj_E0_sz;   //!
   TBranch        *b_detmc_traj_E0;   //!
   TBranch        *b_detmc_traj_Ef_sz;   //!
   TBranch        *b_detmc_traj_Ef;   //!
   TBranch        *b_detmc_traj_preEf_sz;   //!
   TBranch        *b_detmc_traj_preEf;   //!
   TBranch        *b_detmc_traj_prepxf_sz;   //!
   TBranch        *b_detmc_traj_prepxf;   //!
   TBranch        *b_detmc_traj_prepyf_sz;   //!
   TBranch        *b_detmc_traj_prepyf;   //!
   TBranch        *b_detmc_traj_prepzf_sz;   //!
   TBranch        *b_detmc_traj_prepzf;   //!
   TBranch        *b_detmc_traj_px0_sz;   //!
   TBranch        *b_detmc_traj_px0;   //!
   TBranch        *b_detmc_traj_pxf_sz;   //!
   TBranch        *b_detmc_traj_pxf;   //!
   TBranch        *b_detmc_traj_py0_sz;   //!
   TBranch        *b_detmc_traj_py0;   //!
   TBranch        *b_detmc_traj_pyf_sz;   //!
   TBranch        *b_detmc_traj_pyf;   //!
   TBranch        *b_detmc_traj_pz0_sz;   //!
   TBranch        *b_detmc_traj_pz0;   //!
   TBranch        *b_detmc_traj_pzf_sz;   //!
   TBranch        *b_detmc_traj_pzf;   //!
   TBranch        *b_detmc_traj_t0_sz;   //!
   TBranch        *b_detmc_traj_t0;   //!
   TBranch        *b_detmc_traj_tf_sz;   //!
   TBranch        *b_detmc_traj_tf;   //!
   TBranch        *b_detmc_traj_x0_sz;   //!
   TBranch        *b_detmc_traj_x0;   //!
   TBranch        *b_detmc_traj_xf_sz;   //!
   TBranch        *b_detmc_traj_xf;   //!
   TBranch        *b_detmc_traj_y0_sz;   //!
   TBranch        *b_detmc_traj_y0;   //!
   TBranch        *b_detmc_traj_yf_sz;   //!
   TBranch        *b_detmc_traj_yf;   //!
   TBranch        *b_detmc_traj_z0_sz;   //!
   TBranch        *b_detmc_traj_z0;   //!
   TBranch        *b_detmc_traj_zf_sz;   //!
   TBranch        *b_detmc_traj_zf;   //!
   TBranch        *b_recovtx_uvec_sz;   //!
   TBranch        *b_recovtx_uvec;   //!
   TBranch        *b_recovtx_vvec_sz;   //!
   TBranch        *b_recovtx_vvec;   //!
   TBranch        *b_recovtx_xvec_sz;   //!
   TBranch        *b_recovtx_xvec;   //!
   TBranch        *b_recovtx_yvec_sz;   //!
   TBranch        *b_recovtx_yvec;   //!
   TBranch        *b_recovtx_zvec_sz;   //!
   TBranch        *b_recovtx_zvec;   //!
   TBranch        *b_track_cluster_energy_vec_sz;   //!
   TBranch        *b_track_cluster_energy_vec;   //!
   TBranch        *b_track_costheta_vec_sz;   //!
   TBranch        *b_track_costheta_vec;   //!
   TBranch        *b_track_length_vec_sz;   //!
   TBranch        *b_track_length_vec;   //!
   TBranch        *b_track_nodex_vec_sz;   //!
   TBranch        *b_track_nodex_vec;   //!
   TBranch        *b_track_nodey_vec_sz;   //!
   TBranch        *b_track_nodey_vec;   //!
   TBranch        *b_track_nodez_vec_sz;   //!
   TBranch        *b_track_nodez_vec;   //!
   TBranch        *b_track_rev_cluster_energy_vec_sz;   //!
   TBranch        *b_track_rev_cluster_energy_vec;   //!
   TBranch        *b_track_separation_vec_sz;   //!
   TBranch        *b_track_separation_vec;   //!
   TBranch        *b_track_slope_kx_vec_sz;   //!
   TBranch        *b_track_slope_kx_vec;   //!
   TBranch        *b_track_slope_ky_vec_sz;   //!
   TBranch        *b_track_slope_ky_vec;   //!
   TBranch        *b_track_start_uvec_sz;   //!
   TBranch        *b_track_start_uvec;   //!
   TBranch        *b_track_start_vvec_sz;   //!
   TBranch        *b_track_start_vvec;   //!
   TBranch        *b_track_start_xvec_sz;   //!
   TBranch        *b_track_start_xvec;   //!
   TBranch        *b_track_start_yvec_sz;   //!
   TBranch        *b_track_start_yvec;   //!
   TBranch        *b_track_start_zvec_sz;   //!
   TBranch        *b_track_start_zvec;   //!
   TBranch        *b_track_stop_uvec_sz;   //!
   TBranch        *b_track_stop_uvec;   //!
   TBranch        *b_track_stop_vvec_sz;   //!
   TBranch        *b_track_stop_vvec;   //!
   TBranch        *b_track_stop_xvec_sz;   //!
   TBranch        *b_track_stop_xvec;   //!
   TBranch        *b_track_stop_yvec_sz;   //!
   TBranch        *b_track_stop_yvec;   //!
   TBranch        *b_track_stop_zvec_sz;   //!
   TBranch        *b_track_stop_zvec;   //!
   TBranch        *b_track_theta_vec_sz;   //!
   TBranch        *b_track_theta_vec;   //!
   TBranch        *b_track_thetax_vec_sz;   //!
   TBranch        *b_track_thetax_vec;   //!
   TBranch        *b_track_thetay_vec_sz;   //!
   TBranch        *b_track_thetay_vec;   //!
   TBranch        *b_track_truth_fraction_vec1_sz;   //!
   TBranch        *b_track_truth_fraction_vec1;   //!
   TBranch        *b_track_truth_fraction_vec2_sz;   //!
   TBranch        *b_track_truth_fraction_vec2;   //!
   TBranch        *b_track_truth_fraction_vec3_sz;   //!
   TBranch        *b_track_truth_fraction_vec3;   //!
   TBranch        *b_track_truth_shared_vec_sz;   //!
   TBranch        *b_track_truth_shared_vec;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_is_fiducial;   //!
   TBranch        *b_truth_pass_plausible;   //!
   TBranch        *b_truth_is_signal;   //!
   TBranch        *b_truth_nneutron;   //!
   TBranch        *b_truth_nothermeson;   //!
   TBranch        *b_truth_npi0;   //!
   TBranch        *b_truth_npim;   //!
   TBranch        *b_truth_npip;   //!
   TBranch        *b_truth_nproton;   //!
   TBranch        *b_truth_one_pi0;   //!
   TBranch        *b_truth_one_pim;   //!
   TBranch        *b_truth_one_pip;   //!
   TBranch        *b_truth_fslepton_E;   //!
   TBranch        *b_truth_fslepton_P;   //!
   TBranch        *b_truth_fslepton_T;   //!
   TBranch        *b_truth_fslepton_phi;   //!
   TBranch        *b_truth_fslepton_theta;   //!
   TBranch        *b_truth_fslepton_theta_x;   //!
   TBranch        *b_truth_fslepton_theta_y;   //!
   TBranch        *b_dEdxAnaTool_nuFlavor;   //!
   TBranch        *b_dEdxAnaTool_nuHelicity;   //!
   TBranch        *b_dEdxAnaTool_intCurrent;   //!
   TBranch        *b_dEdxAnaTool_intType;   //!
   TBranch        *b_dEdxAnaTool_E;   //!
   TBranch        *b_dEdxAnaTool_Q2;   //!
   TBranch        *b_dEdxAnaTool_x;   //!
   TBranch        *b_dEdxAnaTool_y;   //!
   TBranch        *b_dEdxAnaTool_W;   //!
   TBranch        *b_dEdxAnaTool_score;   //!
   TBranch        *b_dEdxAnaTool_leptonE;   //!
   TBranch        *b_dEdxAnaTool_vtx;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_subrun;   //!
   TBranch        *b_ev_detector;   //!
   TBranch        *b_ev_triggerType;   //!
   TBranch        *b_ev_gate;   //!
   TBranch        *b_ev_global_gate;   //!
   TBranch        *b_ev_gps_time_sec;   //!
   TBranch        *b_ev_gps_time_usec;   //!
   TBranch        *b_mc_run;   //!
   TBranch        *b_mc_subrun;   //!
   TBranch        *b_mc_nInteractions;   //!
   TBranch        *b_mc_MIState;   //!
   TBranch        *b_mc_pot;   //!
   TBranch        *b_mc_beamConfig;   //!
   TBranch        *b_mc_processType;   //!
   TBranch        *b_mc_nthEvtInSpill;   //!
   TBranch        *b_mc_nthEvtInFile;   //!
   TBranch        *b_mc_intType;   //!
   TBranch        *b_mc_current;   //!
   TBranch        *b_mc_charm;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_XSec;   //!
   TBranch        *b_mc_diffXSec;   //!
   TBranch        *b_mc_incoming;   //!
   TBranch        *b_mc_fluxDriverProb;   //!
   TBranch        *b_mc_targetNucleus;   //!
   TBranch        *b_mc_targetZ;   //!
   TBranch        *b_mc_targetA;   //!
   TBranch        *b_mc_targetNucleon;   //!
   TBranch        *b_mc_struckQuark;   //!
   TBranch        *b_mc_seaQuark;   //!
   TBranch        *b_mc_resID;   //!
   TBranch        *b_mc_primaryLepton;   //!
   TBranch        *b_mc_incomingE;   //!
   TBranch        *b_mc_Bjorkenx;   //!
   TBranch        *b_mc_Bjorkeny;   //!
   TBranch        *b_mc_Q2;   //!
   TBranch        *b_mc_nuT;   //!
   TBranch        *b_mc_w;   //!
   TBranch        *b_mc_vtx;   //!
   TBranch        *b_mc_incomingPartVec;   //!
   TBranch        *b_mc_initNucVec;   //!
   TBranch        *b_mc_primFSLepton;   //!
   TBranch        *b_mc_nFSPart;   //!
   TBranch        *b_mc_FSPartPx;   //!
   TBranch        *b_mc_FSPartPy;   //!
   TBranch        *b_mc_FSPartPz;   //!
   TBranch        *b_mc_FSPartE;   //!
   TBranch        *b_mc_FSPartPDG;   //!
   TBranch        *b_mc_er_nPart;   //!
   TBranch        *b_mc_er_ID;   //!
   TBranch        *b_mc_er_status;   //!
   TBranch        *b_mc_er_posInNucX;   //!
   TBranch        *b_mc_er_posInNucY;   //!
   TBranch        *b_mc_er_posInNucZ;   //!
   TBranch        *b_mc_er_Px;   //!
   TBranch        *b_mc_er_Py;   //!
   TBranch        *b_mc_er_Pz;   //!
   TBranch        *b_mc_er_E;   //!
   TBranch        *b_mc_er_FD;   //!
   TBranch        *b_mc_er_LD;   //!
   TBranch        *b_mc_er_mother;   //!
   TBranch        *b_mc_fr_nNuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuParentID;   //!
   TBranch        *b_mc_fr_decMode;   //!
   TBranch        *b_mc_fr_primProtonVtx;   //!
   TBranch        *b_mc_fr_primProtonP;   //!
   TBranch        *b_mc_fr_nuParentDecVtx;   //!
   TBranch        *b_mc_fr_nuParentProdVtx;   //!
   TBranch        *b_mc_fr_nuParentProdP;   //!
   TBranch        *b_mc_cvweight_total;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_mc_cvweight_totalFlux;   //!
   TBranch        *b_mc_cvweight_totalXsec;   //!
   TBranch        *b_mc_ppfx1_cvweight;   //!
   TBranch        *b_mc_hornCurrent_cvweight;   //!
   TBranch        *b_mc_gen1_cvweight_total;   //!
   TBranch        *b_gen1_wgt;   //!
   TBranch        *b_mc_gen1_cvweight_totalFlux;   //!
   TBranch        *b_mc_gen1_cvweight_NA49;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49;   //!
   TBranch        *b_mc_wgt_Norm_sz;   //!
   TBranch        *b_mc_wgt_Norm;   //!
   TBranch        *b_mc_wgt_ppfx1_Total_sz;   //!
   TBranch        *b_mc_wgt_ppfx1_Total;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_prong_nParticles;   //!
   TBranch        *b_prong_part_score;   //!
   TBranch        *b_prong_part_mass;   //!
   TBranch        *b_prong_part_charge;   //!
   TBranch        *b_prong_part_pid;   //!
   TBranch        *b_prong_part_E;   //!
   TBranch        *b_prong_part_pos;   //!

   selector_dedx_calo(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~selector_dedx_calo() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   void Loop();

   ClassDef(selector_dedx_calo,0);
};

#endif

#ifdef selector_dedx_calo_cxx
void selector_dedx_calo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   prong_part_E = 0;
   prong_part_pos = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("physEvtNum", &physEvtNum, &b_physEvtNum);
   fChain->SetBranchAddress("n_hyps", &n_hyps, &b_n_hyps);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("primaryPart", &primaryPart, &b_primaryPart);
   fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
   fChain->SetBranchAddress("slice_numbers", slice_numbers, &b_slice_numbers);
   fChain->SetBranchAddress("shared_slice", &shared_slice, &b_shared_slice);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxErr", vtxErr, &b_vtxErr);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("found_truth", &found_truth, &b_found_truth);
   fChain->SetBranchAddress("phys_front_activity", &phys_front_activity, &b_phys_front_activity);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_is_rockmuon_consistent", &phys_energy_in_road_upstream_is_rockmuon_consistent, &b_phys_energy_in_road_upstream_is_rockmuon_consistent);
   fChain->SetBranchAddress("rock_muons_removed", &rock_muons_removed, &b_rock_muons_removed);
   fChain->SetBranchAddress("minos_track_match", &minos_track_match, &b_minos_track_match);
   fChain->SetBranchAddress("minos_stub_match", &minos_stub_match, &b_minos_stub_match);
   fChain->SetBranchAddress("unknown_helicity", &unknown_helicity, &b_unknown_helicity);
   fChain->SetBranchAddress("minos_track_inside_partial_plane", &minos_track_inside_partial_plane, &b_minos_track_inside_partial_plane);
   fChain->SetBranchAddress("prim_vtx_has_misassigned_track_direction", &prim_vtx_has_misassigned_track_direction, &b_prim_vtx_has_misassigned_track_direction);
   fChain->SetBranchAddress("prim_vtx_has_broken_track", &prim_vtx_has_broken_track, &b_prim_vtx_has_broken_track);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("detmc_ntrajectory", &detmc_ntrajectory, &b_detmc_ntrajectory);
   fChain->SetBranchAddress("detmc_ntrajectory2", &detmc_ntrajectory2, &b_detmc_ntrajectory2);
   fChain->SetBranchAddress("llr_offset_pion", &llr_offset_pion, &b_llr_offset_pion);
   fChain->SetBranchAddress("llr_offset_proton", &llr_offset_proton, &b_llr_offset_proton);
   fChain->SetBranchAddress("minos_trk_end_plane", &minos_trk_end_plane, &b_minos_trk_end_plane);
   fChain->SetBranchAddress("minos_trk_is_contained", &minos_trk_is_contained, &b_minos_trk_is_contained);
   fChain->SetBranchAddress("minos_trk_is_ok", &minos_trk_is_ok, &b_minos_trk_is_ok);
   fChain->SetBranchAddress("minos_trk_quality", &minos_trk_quality, &b_minos_trk_quality);
   fChain->SetBranchAddress("minos_trk_used_curvature", &minos_trk_used_curvature, &b_minos_trk_used_curvature);
   fChain->SetBranchAddress("minos_trk_used_range", &minos_trk_used_range, &b_minos_trk_used_range);
   fChain->SetBranchAddress("ncluster_all", &ncluster_all, &b_ncluster_all);
   fChain->SetBranchAddress("ncluster_used", &ncluster_used, &b_ncluster_used);
   fChain->SetBranchAddress("nprong", &nprong, &b_nprong);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("prong_ntrack", &prong_ntrack, &b_prong_ntrack);
   fChain->SetBranchAddress("recovtx_ntrack", &recovtx_ntrack, &b_recovtx_ntrack);
   fChain->SetBranchAddress("recovtx_ntrack2", &recovtx_ntrack2, &b_recovtx_ntrack2);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("llr_score_pion", &llr_score_pion, &b_llr_score_pion);
   fChain->SetBranchAddress("llr_score_proton", &llr_score_proton, &b_llr_score_proton);
   fChain->SetBranchAddress("minos_trk_eqp", &minos_trk_eqp, &b_minos_trk_eqp);
   fChain->SetBranchAddress("minos_trk_fit_pass", &minos_trk_fit_pass, &b_minos_trk_fit_pass);
   fChain->SetBranchAddress("minos_trk_p", &minos_trk_p, &b_minos_trk_p);
   fChain->SetBranchAddress("minos_trk_p_curvature", &minos_trk_p_curvature, &b_minos_trk_p_curvature);
   fChain->SetBranchAddress("minos_trk_p_range", &minos_trk_p_range, &b_minos_trk_p_range);
   fChain->SetBranchAddress("minos_trk_qp", &minos_trk_qp, &b_minos_trk_qp);
   fChain->SetBranchAddress("muon_E", &muon_E, &b_muon_E);
   fChain->SetBranchAddress("muon_dp", &muon_dp, &b_muon_dp);
   fChain->SetBranchAddress("muon_p", &muon_p, &b_muon_p);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_px", &muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon_py", &muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon_pz", &muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("muon_thetax", &muon_thetax, &b_muon_thetax);
   fChain->SetBranchAddress("muon_thetay", &muon_thetay, &b_muon_thetay);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("pim_gl", &pim_gl, &b_pim_gl);
   fChain->SetBranchAddress("pim_kinetic", &pim_kinetic, &b_pim_kinetic);
   fChain->SetBranchAddress("pim_phi", &pim_phi, &b_pim_phi);
   fChain->SetBranchAddress("pim_pienergy", &pim_pienergy, &b_pim_pienergy);
   fChain->SetBranchAddress("pim_pimom", &pim_pimom, &b_pim_pimom);
   fChain->SetBranchAddress("pim_pl", &pim_pl, &b_pim_pl);
   fChain->SetBranchAddress("pim_pl_area", &pim_pl_area, &b_pim_pl_area);
   fChain->SetBranchAddress("pim_theta", &pim_theta, &b_pim_theta);
   fChain->SetBranchAddress("pim_thetax", &pim_thetax, &b_pim_thetax);
   fChain->SetBranchAddress("pim_thetay", &pim_thetay, &b_pim_thetay);
   fChain->SetBranchAddress("pip_gl", &pip_gl, &b_pip_gl);
   fChain->SetBranchAddress("pip_kinetic", &pip_kinetic, &b_pip_kinetic);
   fChain->SetBranchAddress("pip_phi", &pip_phi, &b_pip_phi);
   fChain->SetBranchAddress("pip_pienergy", &pip_pienergy, &b_pip_pienergy);
   fChain->SetBranchAddress("pip_pimom", &pip_pimom, &b_pip_pimom);
   fChain->SetBranchAddress("pip_pl", &pip_pl, &b_pip_pl);
   fChain->SetBranchAddress("pip_pl_area", &pip_pl_area, &b_pip_pl_area);
   fChain->SetBranchAddress("pip_theta", &pip_theta, &b_pip_theta);
   fChain->SetBranchAddress("pip_thetax", &pip_thetax, &b_pip_thetax);
   fChain->SetBranchAddress("pip_thetay", &pip_thetay, &b_pip_thetay);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("prong_hypo_pion_E", &prong_hypo_pion_E, &b_prong_hypo_pion_E);
   fChain->SetBranchAddress("prong_hypo_pion_charge", &prong_hypo_pion_charge, &b_prong_hypo_pion_charge);
   fChain->SetBranchAddress("prong_hypo_pion_mass", &prong_hypo_pion_mass, &b_prong_hypo_pion_mass);
   fChain->SetBranchAddress("prong_hypo_pion_p", &prong_hypo_pion_p, &b_prong_hypo_pion_p);
   fChain->SetBranchAddress("prong_hypo_pion_px", &prong_hypo_pion_px, &b_prong_hypo_pion_px);
   fChain->SetBranchAddress("prong_hypo_pion_py", &prong_hypo_pion_py, &b_prong_hypo_pion_py);
   fChain->SetBranchAddress("prong_hypo_pion_pz", &prong_hypo_pion_pz, &b_prong_hypo_pion_pz);
   fChain->SetBranchAddress("prong_hypo_pion_score", &prong_hypo_pion_score, &b_prong_hypo_pion_score);
   fChain->SetBranchAddress("prong_hypo_proton_E", &prong_hypo_proton_E, &b_prong_hypo_proton_E);
   fChain->SetBranchAddress("prong_hypo_proton_charge", &prong_hypo_proton_charge, &b_prong_hypo_proton_charge);
   fChain->SetBranchAddress("prong_hypo_proton_mass", &prong_hypo_proton_mass, &b_prong_hypo_proton_mass);
   fChain->SetBranchAddress("prong_hypo_proton_p", &prong_hypo_proton_p, &b_prong_hypo_proton_p);
   fChain->SetBranchAddress("prong_hypo_proton_px", &prong_hypo_proton_px, &b_prong_hypo_proton_px);
   fChain->SetBranchAddress("prong_hypo_proton_py", &prong_hypo_proton_py, &b_prong_hypo_proton_py);
   fChain->SetBranchAddress("prong_hypo_proton_pz", &prong_hypo_proton_pz, &b_prong_hypo_proton_pz);
   fChain->SetBranchAddress("prong_hypo_proton_score", &prong_hypo_proton_score, &b_prong_hypo_proton_score);
   fChain->SetBranchAddress("prong_length", &prong_length, &b_prong_length);
   fChain->SetBranchAddress("prong_length_area", &prong_length_area, &b_prong_length_area);
   fChain->SetBranchAddress("prong_start_u", &prong_start_u, &b_prong_start_u);
   fChain->SetBranchAddress("prong_start_v", &prong_start_v, &b_prong_start_v);
   fChain->SetBranchAddress("prong_start_x", &prong_start_x, &b_prong_start_x);
   fChain->SetBranchAddress("prong_start_y", &prong_start_y, &b_prong_start_y);
   fChain->SetBranchAddress("prong_start_z", &prong_start_z, &b_prong_start_z);
   fChain->SetBranchAddress("prong_stop_u", &prong_stop_u, &b_prong_stop_u);
   fChain->SetBranchAddress("prong_stop_v", &prong_stop_v, &b_prong_stop_v);
   fChain->SetBranchAddress("prong_stop_x", &prong_stop_x, &b_prong_stop_x);
   fChain->SetBranchAddress("prong_stop_y", &prong_stop_y, &b_prong_stop_y);
   fChain->SetBranchAddress("prong_stop_z", &prong_stop_z, &b_prong_stop_z);
   fChain->SetBranchAddress("cluster_used_vec_sz", &cluster_used_vec_sz, &b_cluster_used_vec_sz);
   fChain->SetBranchAddress("cluster_used_vec", cluster_used_vec, &b_cluster_used_vec);
   fChain->SetBranchAddress("cluster_view_vec_sz", &cluster_view_vec_sz, &b_cluster_view_vec_sz);
   fChain->SetBranchAddress("cluster_view_vec", cluster_view_vec, &b_cluster_view_vec);
   fChain->SetBranchAddress("detmc_traj_id_sz", &detmc_traj_id_sz, &b_detmc_traj_id_sz);
   fChain->SetBranchAddress("detmc_traj_id", detmc_traj_id, &b_detmc_traj_id);
   fChain->SetBranchAddress("detmc_traj_mother_sz", &detmc_traj_mother_sz, &b_detmc_traj_mother_sz);
   fChain->SetBranchAddress("detmc_traj_mother", detmc_traj_mother, &b_detmc_traj_mother);
   fChain->SetBranchAddress("detmc_traj_pdg_sz", &detmc_traj_pdg_sz, &b_detmc_traj_pdg_sz);
   fChain->SetBranchAddress("detmc_traj_pdg", detmc_traj_pdg, &b_detmc_traj_pdg);
   fChain->SetBranchAddress("detmc_traj_proc_sz", &detmc_traj_proc_sz, &b_detmc_traj_proc_sz);
   fChain->SetBranchAddress("detmc_traj_proc", detmc_traj_proc, &b_detmc_traj_proc);
   fChain->SetBranchAddress("detmc_traj_status_sz", &detmc_traj_status_sz, &b_detmc_traj_status_sz);
   fChain->SetBranchAddress("detmc_traj_status", detmc_traj_status, &b_detmc_traj_status);
   fChain->SetBranchAddress("track_module_vec_sz", &track_module_vec_sz, &b_track_module_vec_sz);
   fChain->SetBranchAddress("track_module_vec", track_module_vec, &b_track_module_vec);
   fChain->SetBranchAddress("track_ndigit_vec_sz", &track_ndigit_vec_sz, &b_track_ndigit_vec_sz);
   fChain->SetBranchAddress("track_ndigit_vec", track_ndigit_vec, &b_track_ndigit_vec);
   fChain->SetBranchAddress("track_node_vec_sz", &track_node_vec_sz, &b_track_node_vec_sz);
   fChain->SetBranchAddress("track_node_vec", track_node_vec, &b_track_node_vec);
   fChain->SetBranchAddress("track_plane_vec_sz", &track_plane_vec_sz, &b_track_plane_vec_sz);
   fChain->SetBranchAddress("track_plane_vec", track_plane_vec, &b_track_plane_vec);
   fChain->SetBranchAddress("track_truth_count_vec_sz", &track_truth_count_vec_sz, &b_track_truth_count_vec_sz);
   fChain->SetBranchAddress("track_truth_count_vec", track_truth_count_vec, &b_track_truth_count_vec);
   fChain->SetBranchAddress("track_truth_pdg_vec1_sz", &track_truth_pdg_vec1_sz, &b_track_truth_pdg_vec1_sz);
   fChain->SetBranchAddress("track_truth_pdg_vec1", track_truth_pdg_vec1, &b_track_truth_pdg_vec1);
   fChain->SetBranchAddress("track_truth_pdg_vec2_sz", &track_truth_pdg_vec2_sz, &b_track_truth_pdg_vec2_sz);
   fChain->SetBranchAddress("track_truth_pdg_vec2", track_truth_pdg_vec2, &b_track_truth_pdg_vec2);
   fChain->SetBranchAddress("track_truth_pdg_vec3_sz", &track_truth_pdg_vec3_sz, &b_track_truth_pdg_vec3_sz);
   fChain->SetBranchAddress("track_truth_pdg_vec3", track_truth_pdg_vec3, &b_track_truth_pdg_vec3);
   fChain->SetBranchAddress("cluster_energy_vec_sz", &cluster_energy_vec_sz, &b_cluster_energy_vec_sz);
   fChain->SetBranchAddress("cluster_energy_vec", cluster_energy_vec, &b_cluster_energy_vec);
   fChain->SetBranchAddress("cluster_time_vec_sz", &cluster_time_vec_sz, &b_cluster_time_vec_sz);
   fChain->SetBranchAddress("cluster_time_vec", cluster_time_vec, &b_cluster_time_vec);
   fChain->SetBranchAddress("cluster_x_vec_sz", &cluster_x_vec_sz, &b_cluster_x_vec_sz);
   fChain->SetBranchAddress("cluster_x_vec", cluster_x_vec, &b_cluster_x_vec);
   fChain->SetBranchAddress("cluster_z_vec_sz", &cluster_z_vec_sz, &b_cluster_z_vec_sz);
   fChain->SetBranchAddress("cluster_z_vec", cluster_z_vec, &b_cluster_z_vec);
   fChain->SetBranchAddress("detmc_traj_E0_sz", &detmc_traj_E0_sz, &b_detmc_traj_E0_sz);
   fChain->SetBranchAddress("detmc_traj_E0", detmc_traj_E0, &b_detmc_traj_E0);
   fChain->SetBranchAddress("detmc_traj_Ef_sz", &detmc_traj_Ef_sz, &b_detmc_traj_Ef_sz);
   fChain->SetBranchAddress("detmc_traj_Ef", detmc_traj_Ef, &b_detmc_traj_Ef);
   fChain->SetBranchAddress("detmc_traj_preEf_sz", &detmc_traj_preEf_sz, &b_detmc_traj_preEf_sz);
   fChain->SetBranchAddress("detmc_traj_preEf", detmc_traj_preEf, &b_detmc_traj_preEf);
   fChain->SetBranchAddress("detmc_traj_prepxf_sz", &detmc_traj_prepxf_sz, &b_detmc_traj_prepxf_sz);
   fChain->SetBranchAddress("detmc_traj_prepxf", detmc_traj_prepxf, &b_detmc_traj_prepxf);
   fChain->SetBranchAddress("detmc_traj_prepyf_sz", &detmc_traj_prepyf_sz, &b_detmc_traj_prepyf_sz);
   fChain->SetBranchAddress("detmc_traj_prepyf", detmc_traj_prepyf, &b_detmc_traj_prepyf);
   fChain->SetBranchAddress("detmc_traj_prepzf_sz", &detmc_traj_prepzf_sz, &b_detmc_traj_prepzf_sz);
   fChain->SetBranchAddress("detmc_traj_prepzf", detmc_traj_prepzf, &b_detmc_traj_prepzf);
   fChain->SetBranchAddress("detmc_traj_px0_sz", &detmc_traj_px0_sz, &b_detmc_traj_px0_sz);
   fChain->SetBranchAddress("detmc_traj_px0", detmc_traj_px0, &b_detmc_traj_px0);
   fChain->SetBranchAddress("detmc_traj_pxf_sz", &detmc_traj_pxf_sz, &b_detmc_traj_pxf_sz);
   fChain->SetBranchAddress("detmc_traj_pxf", detmc_traj_pxf, &b_detmc_traj_pxf);
   fChain->SetBranchAddress("detmc_traj_py0_sz", &detmc_traj_py0_sz, &b_detmc_traj_py0_sz);
   fChain->SetBranchAddress("detmc_traj_py0", detmc_traj_py0, &b_detmc_traj_py0);
   fChain->SetBranchAddress("detmc_traj_pyf_sz", &detmc_traj_pyf_sz, &b_detmc_traj_pyf_sz);
   fChain->SetBranchAddress("detmc_traj_pyf", detmc_traj_pyf, &b_detmc_traj_pyf);
   fChain->SetBranchAddress("detmc_traj_pz0_sz", &detmc_traj_pz0_sz, &b_detmc_traj_pz0_sz);
   fChain->SetBranchAddress("detmc_traj_pz0", detmc_traj_pz0, &b_detmc_traj_pz0);
   fChain->SetBranchAddress("detmc_traj_pzf_sz", &detmc_traj_pzf_sz, &b_detmc_traj_pzf_sz);
   fChain->SetBranchAddress("detmc_traj_pzf", detmc_traj_pzf, &b_detmc_traj_pzf);
   fChain->SetBranchAddress("detmc_traj_t0_sz", &detmc_traj_t0_sz, &b_detmc_traj_t0_sz);
   fChain->SetBranchAddress("detmc_traj_t0", detmc_traj_t0, &b_detmc_traj_t0);
   fChain->SetBranchAddress("detmc_traj_tf_sz", &detmc_traj_tf_sz, &b_detmc_traj_tf_sz);
   fChain->SetBranchAddress("detmc_traj_tf", detmc_traj_tf, &b_detmc_traj_tf);
   fChain->SetBranchAddress("detmc_traj_x0_sz", &detmc_traj_x0_sz, &b_detmc_traj_x0_sz);
   fChain->SetBranchAddress("detmc_traj_x0", detmc_traj_x0, &b_detmc_traj_x0);
   fChain->SetBranchAddress("detmc_traj_xf_sz", &detmc_traj_xf_sz, &b_detmc_traj_xf_sz);
   fChain->SetBranchAddress("detmc_traj_xf", detmc_traj_xf, &b_detmc_traj_xf);
   fChain->SetBranchAddress("detmc_traj_y0_sz", &detmc_traj_y0_sz, &b_detmc_traj_y0_sz);
   fChain->SetBranchAddress("detmc_traj_y0", detmc_traj_y0, &b_detmc_traj_y0);
   fChain->SetBranchAddress("detmc_traj_yf_sz", &detmc_traj_yf_sz, &b_detmc_traj_yf_sz);
   fChain->SetBranchAddress("detmc_traj_yf", detmc_traj_yf, &b_detmc_traj_yf);
   fChain->SetBranchAddress("detmc_traj_z0_sz", &detmc_traj_z0_sz, &b_detmc_traj_z0_sz);
   fChain->SetBranchAddress("detmc_traj_z0", detmc_traj_z0, &b_detmc_traj_z0);
   fChain->SetBranchAddress("detmc_traj_zf_sz", &detmc_traj_zf_sz, &b_detmc_traj_zf_sz);
   fChain->SetBranchAddress("detmc_traj_zf", detmc_traj_zf, &b_detmc_traj_zf);
   fChain->SetBranchAddress("recovtx_uvec_sz", &recovtx_uvec_sz, &b_recovtx_uvec_sz);
   fChain->SetBranchAddress("recovtx_uvec", recovtx_uvec, &b_recovtx_uvec);
   fChain->SetBranchAddress("recovtx_vvec_sz", &recovtx_vvec_sz, &b_recovtx_vvec_sz);
   fChain->SetBranchAddress("recovtx_vvec", recovtx_vvec, &b_recovtx_vvec);
   fChain->SetBranchAddress("recovtx_xvec_sz", &recovtx_xvec_sz, &b_recovtx_xvec_sz);
   fChain->SetBranchAddress("recovtx_xvec", recovtx_xvec, &b_recovtx_xvec);
   fChain->SetBranchAddress("recovtx_yvec_sz", &recovtx_yvec_sz, &b_recovtx_yvec_sz);
   fChain->SetBranchAddress("recovtx_yvec", recovtx_yvec, &b_recovtx_yvec);
   fChain->SetBranchAddress("recovtx_zvec_sz", &recovtx_zvec_sz, &b_recovtx_zvec_sz);
   fChain->SetBranchAddress("recovtx_zvec", recovtx_zvec, &b_recovtx_zvec);
   fChain->SetBranchAddress("track_cluster_energy_vec_sz", &track_cluster_energy_vec_sz, &b_track_cluster_energy_vec_sz);
   fChain->SetBranchAddress("track_cluster_energy_vec", track_cluster_energy_vec, &b_track_cluster_energy_vec);
   fChain->SetBranchAddress("track_costheta_vec_sz", &track_costheta_vec_sz, &b_track_costheta_vec_sz);
   fChain->SetBranchAddress("track_costheta_vec", track_costheta_vec, &b_track_costheta_vec);
   fChain->SetBranchAddress("track_length_vec_sz", &track_length_vec_sz, &b_track_length_vec_sz);
   fChain->SetBranchAddress("track_length_vec", track_length_vec, &b_track_length_vec);
   fChain->SetBranchAddress("track_nodex_vec_sz", &track_nodex_vec_sz, &b_track_nodex_vec_sz);
   fChain->SetBranchAddress("track_nodex_vec", track_nodex_vec, &b_track_nodex_vec);
   fChain->SetBranchAddress("track_nodey_vec_sz", &track_nodey_vec_sz, &b_track_nodey_vec_sz);
   fChain->SetBranchAddress("track_nodey_vec", track_nodey_vec, &b_track_nodey_vec);
   fChain->SetBranchAddress("track_nodez_vec_sz", &track_nodez_vec_sz, &b_track_nodez_vec_sz);
   fChain->SetBranchAddress("track_nodez_vec", track_nodez_vec, &b_track_nodez_vec);
   fChain->SetBranchAddress("track_rev_cluster_energy_vec_sz", &track_rev_cluster_energy_vec_sz, &b_track_rev_cluster_energy_vec_sz);
   fChain->SetBranchAddress("track_rev_cluster_energy_vec", track_rev_cluster_energy_vec, &b_track_rev_cluster_energy_vec);
   fChain->SetBranchAddress("track_separation_vec_sz", &track_separation_vec_sz, &b_track_separation_vec_sz);
   fChain->SetBranchAddress("track_separation_vec", track_separation_vec, &b_track_separation_vec);
   fChain->SetBranchAddress("track_slope_kx_vec_sz", &track_slope_kx_vec_sz, &b_track_slope_kx_vec_sz);
   fChain->SetBranchAddress("track_slope_kx_vec", track_slope_kx_vec, &b_track_slope_kx_vec);
   fChain->SetBranchAddress("track_slope_ky_vec_sz", &track_slope_ky_vec_sz, &b_track_slope_ky_vec_sz);
   fChain->SetBranchAddress("track_slope_ky_vec", track_slope_ky_vec, &b_track_slope_ky_vec);
   fChain->SetBranchAddress("track_start_uvec_sz", &track_start_uvec_sz, &b_track_start_uvec_sz);
   fChain->SetBranchAddress("track_start_uvec", track_start_uvec, &b_track_start_uvec);
   fChain->SetBranchAddress("track_start_vvec_sz", &track_start_vvec_sz, &b_track_start_vvec_sz);
   fChain->SetBranchAddress("track_start_vvec", track_start_vvec, &b_track_start_vvec);
   fChain->SetBranchAddress("track_start_xvec_sz", &track_start_xvec_sz, &b_track_start_xvec_sz);
   fChain->SetBranchAddress("track_start_xvec", track_start_xvec, &b_track_start_xvec);
   fChain->SetBranchAddress("track_start_yvec_sz", &track_start_yvec_sz, &b_track_start_yvec_sz);
   fChain->SetBranchAddress("track_start_yvec", track_start_yvec, &b_track_start_yvec);
   fChain->SetBranchAddress("track_start_zvec_sz", &track_start_zvec_sz, &b_track_start_zvec_sz);
   fChain->SetBranchAddress("track_start_zvec", track_start_zvec, &b_track_start_zvec);
   fChain->SetBranchAddress("track_stop_uvec_sz", &track_stop_uvec_sz, &b_track_stop_uvec_sz);
   fChain->SetBranchAddress("track_stop_uvec", track_stop_uvec, &b_track_stop_uvec);
   fChain->SetBranchAddress("track_stop_vvec_sz", &track_stop_vvec_sz, &b_track_stop_vvec_sz);
   fChain->SetBranchAddress("track_stop_vvec", track_stop_vvec, &b_track_stop_vvec);
   fChain->SetBranchAddress("track_stop_xvec_sz", &track_stop_xvec_sz, &b_track_stop_xvec_sz);
   fChain->SetBranchAddress("track_stop_xvec", track_stop_xvec, &b_track_stop_xvec);
   fChain->SetBranchAddress("track_stop_yvec_sz", &track_stop_yvec_sz, &b_track_stop_yvec_sz);
   fChain->SetBranchAddress("track_stop_yvec", track_stop_yvec, &b_track_stop_yvec);
   fChain->SetBranchAddress("track_stop_zvec_sz", &track_stop_zvec_sz, &b_track_stop_zvec_sz);
   fChain->SetBranchAddress("track_stop_zvec", track_stop_zvec, &b_track_stop_zvec);
   fChain->SetBranchAddress("track_theta_vec_sz", &track_theta_vec_sz, &b_track_theta_vec_sz);
   fChain->SetBranchAddress("track_theta_vec", track_theta_vec, &b_track_theta_vec);
   fChain->SetBranchAddress("track_thetax_vec_sz", &track_thetax_vec_sz, &b_track_thetax_vec_sz);
   fChain->SetBranchAddress("track_thetax_vec", track_thetax_vec, &b_track_thetax_vec);
   fChain->SetBranchAddress("track_thetay_vec_sz", &track_thetay_vec_sz, &b_track_thetay_vec_sz);
   fChain->SetBranchAddress("track_thetay_vec", track_thetay_vec, &b_track_thetay_vec);
   fChain->SetBranchAddress("track_truth_fraction_vec1_sz", &track_truth_fraction_vec1_sz, &b_track_truth_fraction_vec1_sz);
   fChain->SetBranchAddress("track_truth_fraction_vec1", track_truth_fraction_vec1, &b_track_truth_fraction_vec1);
   fChain->SetBranchAddress("track_truth_fraction_vec2_sz", &track_truth_fraction_vec2_sz, &b_track_truth_fraction_vec2_sz);
   fChain->SetBranchAddress("track_truth_fraction_vec2", track_truth_fraction_vec2, &b_track_truth_fraction_vec2);
   fChain->SetBranchAddress("track_truth_fraction_vec3_sz", &track_truth_fraction_vec3_sz, &b_track_truth_fraction_vec3_sz);
   fChain->SetBranchAddress("track_truth_fraction_vec3", track_truth_fraction_vec3, &b_track_truth_fraction_vec3);
   fChain->SetBranchAddress("track_truth_shared_vec_sz", &track_truth_shared_vec_sz, &b_track_truth_shared_vec_sz);
   fChain->SetBranchAddress("track_truth_shared_vec", track_truth_shared_vec, &b_track_truth_shared_vec);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_is_fiducial", &truth_is_fiducial, &b_truth_is_fiducial);
   fChain->SetBranchAddress("truth_pass_plausible", &truth_pass_plausible, &b_truth_pass_plausible);
   fChain->SetBranchAddress("truth_is_signal", &truth_is_signal, &b_truth_is_signal);
   fChain->SetBranchAddress("truth_nneutron", &truth_nneutron, &b_truth_nneutron);
   fChain->SetBranchAddress("truth_nothermeson", &truth_nothermeson, &b_truth_nothermeson);
   fChain->SetBranchAddress("truth_npi0", &truth_npi0, &b_truth_npi0);
   fChain->SetBranchAddress("truth_npim", &truth_npim, &b_truth_npim);
   fChain->SetBranchAddress("truth_npip", &truth_npip, &b_truth_npip);
   fChain->SetBranchAddress("truth_nproton", &truth_nproton, &b_truth_nproton);
   fChain->SetBranchAddress("truth_one_pi0", &truth_one_pi0, &b_truth_one_pi0);
   fChain->SetBranchAddress("truth_one_pim", &truth_one_pim, &b_truth_one_pim);
   fChain->SetBranchAddress("truth_one_pip", &truth_one_pip, &b_truth_one_pip);
   fChain->SetBranchAddress("truth_fslepton_E", &truth_fslepton_E, &b_truth_fslepton_E);
   fChain->SetBranchAddress("truth_fslepton_P", &truth_fslepton_P, &b_truth_fslepton_P);
   fChain->SetBranchAddress("truth_fslepton_T", &truth_fslepton_T, &b_truth_fslepton_T);
   fChain->SetBranchAddress("truth_fslepton_phi", &truth_fslepton_phi, &b_truth_fslepton_phi);
   fChain->SetBranchAddress("truth_fslepton_theta", &truth_fslepton_theta, &b_truth_fslepton_theta);
   fChain->SetBranchAddress("truth_fslepton_theta_x", &truth_fslepton_theta_x, &b_truth_fslepton_theta_x);
   fChain->SetBranchAddress("truth_fslepton_theta_y", &truth_fslepton_theta_y, &b_truth_fslepton_theta_y);
   fChain->SetBranchAddress("dEdxAnaTool_nuFlavor", &dEdxAnaTool_nuFlavor, &b_dEdxAnaTool_nuFlavor);
   fChain->SetBranchAddress("dEdxAnaTool_nuHelicity", &dEdxAnaTool_nuHelicity, &b_dEdxAnaTool_nuHelicity);
   fChain->SetBranchAddress("dEdxAnaTool_intCurrent", &dEdxAnaTool_intCurrent, &b_dEdxAnaTool_intCurrent);
   fChain->SetBranchAddress("dEdxAnaTool_intType", &dEdxAnaTool_intType, &b_dEdxAnaTool_intType);
   fChain->SetBranchAddress("dEdxAnaTool_E", &dEdxAnaTool_E, &b_dEdxAnaTool_E);
   fChain->SetBranchAddress("dEdxAnaTool_Q2", &dEdxAnaTool_Q2, &b_dEdxAnaTool_Q2);
   fChain->SetBranchAddress("dEdxAnaTool_x", &dEdxAnaTool_x, &b_dEdxAnaTool_x);
   fChain->SetBranchAddress("dEdxAnaTool_y", &dEdxAnaTool_y, &b_dEdxAnaTool_y);
   fChain->SetBranchAddress("dEdxAnaTool_W", &dEdxAnaTool_W, &b_dEdxAnaTool_W);
   fChain->SetBranchAddress("dEdxAnaTool_score", &dEdxAnaTool_score, &b_dEdxAnaTool_score);
   fChain->SetBranchAddress("dEdxAnaTool_leptonE", dEdxAnaTool_leptonE, &b_dEdxAnaTool_leptonE);
   fChain->SetBranchAddress("dEdxAnaTool_vtx", dEdxAnaTool_vtx, &b_dEdxAnaTool_vtx);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_subrun", &ev_subrun, &b_ev_subrun);
   fChain->SetBranchAddress("ev_detector", &ev_detector, &b_ev_detector);
   fChain->SetBranchAddress("ev_triggerType", &ev_triggerType, &b_ev_triggerType);
   fChain->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
   fChain->SetBranchAddress("ev_global_gate", &ev_global_gate, &b_ev_global_gate);
   fChain->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec, &b_ev_gps_time_sec);
   fChain->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec, &b_ev_gps_time_usec);
   fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
   fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
   fChain->SetBranchAddress("mc_nInteractions", &mc_nInteractions, &b_mc_nInteractions);
   fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
   fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
   fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
   fChain->SetBranchAddress("mc_processType", &mc_processType, &b_mc_processType);
   fChain->SetBranchAddress("mc_nthEvtInSpill", &mc_nthEvtInSpill, &b_mc_nthEvtInSpill);
   fChain->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
   fChain->SetBranchAddress("mc_intType", &mc_intType, &b_mc_intType);
   fChain->SetBranchAddress("mc_current", &mc_current, &b_mc_current);
   fChain->SetBranchAddress("mc_charm", &mc_charm, &b_mc_charm);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("mc_XSec", &mc_XSec, &b_mc_XSec);
   fChain->SetBranchAddress("mc_diffXSec", &mc_diffXSec, &b_mc_diffXSec);
   fChain->SetBranchAddress("mc_incoming", &mc_incoming, &b_mc_incoming);
   fChain->SetBranchAddress("mc_fluxDriverProb", &mc_fluxDriverProb, &b_mc_fluxDriverProb);
   fChain->SetBranchAddress("mc_targetNucleus", &mc_targetNucleus, &b_mc_targetNucleus);
   fChain->SetBranchAddress("mc_targetZ", &mc_targetZ, &b_mc_targetZ);
   fChain->SetBranchAddress("mc_targetA", &mc_targetA, &b_mc_targetA);
   fChain->SetBranchAddress("mc_targetNucleon", &mc_targetNucleon, &b_mc_targetNucleon);
   fChain->SetBranchAddress("mc_struckQuark", &mc_struckQuark, &b_mc_struckQuark);
   fChain->SetBranchAddress("mc_seaQuark", &mc_seaQuark, &b_mc_seaQuark);
   fChain->SetBranchAddress("mc_resID", &mc_resID, &b_mc_resID);
   fChain->SetBranchAddress("mc_primaryLepton", &mc_primaryLepton, &b_mc_primaryLepton);
   fChain->SetBranchAddress("mc_incomingE", &mc_incomingE, &b_mc_incomingE);
   fChain->SetBranchAddress("mc_Bjorkenx", &mc_Bjorkenx, &b_mc_Bjorkenx);
   fChain->SetBranchAddress("mc_Bjorkeny", &mc_Bjorkeny, &b_mc_Bjorkeny);
   fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
   fChain->SetBranchAddress("mc_nuT", &mc_nuT, &b_mc_nuT);
   fChain->SetBranchAddress("mc_w", &mc_w, &b_mc_w);
   fChain->SetBranchAddress("mc_vtx", mc_vtx, &b_mc_vtx);
   fChain->SetBranchAddress("mc_incomingPartVec", mc_incomingPartVec, &b_mc_incomingPartVec);
   fChain->SetBranchAddress("mc_initNucVec", mc_initNucVec, &b_mc_initNucVec);
   fChain->SetBranchAddress("mc_primFSLepton", mc_primFSLepton, &b_mc_primFSLepton);
   fChain->SetBranchAddress("mc_nFSPart", &mc_nFSPart, &b_mc_nFSPart);
   fChain->SetBranchAddress("mc_FSPartPx", mc_FSPartPx, &b_mc_FSPartPx);
   fChain->SetBranchAddress("mc_FSPartPy", mc_FSPartPy, &b_mc_FSPartPy);
   fChain->SetBranchAddress("mc_FSPartPz", mc_FSPartPz, &b_mc_FSPartPz);
   fChain->SetBranchAddress("mc_FSPartE", mc_FSPartE, &b_mc_FSPartE);
   fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG);
   fChain->SetBranchAddress("mc_er_nPart", &mc_er_nPart, &b_mc_er_nPart);
   fChain->SetBranchAddress("mc_er_ID", &mc_er_ID, &b_mc_er_ID);
   fChain->SetBranchAddress("mc_er_status", &mc_er_status, &b_mc_er_status);
   fChain->SetBranchAddress("mc_er_posInNucX", &mc_er_posInNucX, &b_mc_er_posInNucX);
   fChain->SetBranchAddress("mc_er_posInNucY", &mc_er_posInNucY, &b_mc_er_posInNucY);
   fChain->SetBranchAddress("mc_er_posInNucZ", &mc_er_posInNucZ, &b_mc_er_posInNucZ);
   fChain->SetBranchAddress("mc_er_Px", &mc_er_Px, &b_mc_er_Px);
   fChain->SetBranchAddress("mc_er_Py", &mc_er_Py, &b_mc_er_Py);
   fChain->SetBranchAddress("mc_er_Pz", &mc_er_Pz, &b_mc_er_Pz);
   fChain->SetBranchAddress("mc_er_E", &mc_er_E, &b_mc_er_E);
   fChain->SetBranchAddress("mc_er_FD", &mc_er_FD, &b_mc_er_FD);
   fChain->SetBranchAddress("mc_er_LD", &mc_er_LD, &b_mc_er_LD);
   fChain->SetBranchAddress("mc_er_mother", &mc_er_mother, &b_mc_er_mother);
   fChain->SetBranchAddress("mc_fr_nNuAncestorIDs", &mc_fr_nNuAncestorIDs, &b_mc_fr_nNuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuAncestorIDs", &mc_fr_nuAncestorIDs, &b_mc_fr_nuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuParentID", &mc_fr_nuParentID, &b_mc_fr_nuParentID);
   fChain->SetBranchAddress("mc_fr_decMode", &mc_fr_decMode, &b_mc_fr_decMode);
   fChain->SetBranchAddress("mc_fr_primProtonVtx", mc_fr_primProtonVtx, &b_mc_fr_primProtonVtx);
   fChain->SetBranchAddress("mc_fr_primProtonP", mc_fr_primProtonP, &b_mc_fr_primProtonP);
   fChain->SetBranchAddress("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, &b_mc_fr_nuParentDecVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, &b_mc_fr_nuParentProdVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdP", mc_fr_nuParentProdP, &b_mc_fr_nuParentProdP);
   fChain->SetBranchAddress("mc_cvweight_total", &mc_cvweight_total, &b_mc_cvweight_total);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, &b_mc_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, &b_mc_cvweight_totalXsec);
   fChain->SetBranchAddress("mc_ppfx1_cvweight", &mc_ppfx1_cvweight, &b_mc_ppfx1_cvweight);
   fChain->SetBranchAddress("mc_hornCurrent_cvweight", &mc_hornCurrent_cvweight, &b_mc_hornCurrent_cvweight);
   fChain->SetBranchAddress("mc_gen1_cvweight_total", &mc_gen1_cvweight_total, &b_mc_gen1_cvweight_total);
   fChain->SetBranchAddress("gen1_wgt", &gen1_wgt, &b_gen1_wgt);
   fChain->SetBranchAddress("mc_gen1_cvweight_totalFlux", &mc_gen1_cvweight_totalFlux, &b_mc_gen1_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_gen1_cvweight_NA49", &mc_gen1_cvweight_NA49, &b_mc_gen1_cvweight_NA49);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", &mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary_sz", &mc_wgt_gen1_Flux_Tertiary_sz, &b_mc_wgt_gen1_Flux_Tertiary_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary", &mc_wgt_gen1_Flux_Tertiary, &b_mc_wgt_gen1_Flux_Tertiary);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49_sz", &mc_wgt_gen1_Flux_NA49_sz, &b_mc_wgt_gen1_Flux_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49", &mc_wgt_gen1_Flux_NA49, &b_mc_wgt_gen1_Flux_NA49);
   fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
   fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total_sz", &mc_wgt_ppfx1_Total_sz, &b_mc_wgt_ppfx1_Total_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total", &mc_wgt_ppfx1_Total, &b_mc_wgt_ppfx1_Total);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
   fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
   fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
   fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
   fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
//   fChain->SetBranchAddress("prong_part_E", &prong_part_E, &b_prong_part_E);
//   fChain->SetBranchAddress("prong_part_pos", &prong_part_pos, &b_prong_part_pos);
}

Bool_t selector_dedx_calo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef selector_dedx_calo_cxx

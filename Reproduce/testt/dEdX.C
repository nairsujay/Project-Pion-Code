void dEdX()
{
	TH2D* myhist = new TH2D("myhist","Histogram of dE/dx for 6 steps",6,0,6,10,0,40);
	TFile::Open("OutputMINERvA_tagmc_minerva_plastic_target_r4871_Neutrino2016poster_331127_6463088_ProtonCorr0_MuonCorr0_doCut0.root");
	TTree* tree=(TTree*)gDirectory->Get("tree");
	Double_t E0,E1,E2,E3,E4,E5;

	tree->SetBranchAddress("fProtonRecE0",&E0);

	tree->SetBranchAddress("fProtonRecE1",&E1);
	
	tree->SetBranchAddress("fProtonRecE2",&E2);

	tree->SetBranchAddress("fProtonRecE3",&E3);

	tree->SetBranchAddress("fProtonRecE4",&E4);

	tree->SetBranchAddress("fProtonRecE5",&E5);
	for(Int_t ii=0, N=tree->GetEntries(); ii<N;++ii)
	{
		tree->GetEntry(ii);
		myhist->Fill(0,E0);
		myhist->Fill(1,E1);
		myhist->Fill(2,E2);
		myhist->Fill(3,E3);
		myhist->Fill(4,E4);
		myhist->Fill(5,E5);

	}
	myhist->SetStats(kFALSE);
	myhist->SetXTitle("Step");
	myhist->SetYTitle("Energy Deposited/Step");
	myhist->Draw("COLZ");
	c1->Print("/home/jesu2900/project/2017Summer/MINERvA/testt/dEdX.png");

}

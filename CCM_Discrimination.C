void CCM_Discrimination(Float_t CCM_thresh = 0.697)
{
	gROOT->Reset();
	TH1F* h_Spectrum_gamma = new TH1F("h_Spectrum_gamma","SGB 2inch LaBr3 gamma;Energy/keV;Counts",300,1700,2800);
	TH1F* h_Spectrum_alpha = new TH1F("h_Spectrum_alpha","SGB 2inch LaBr3 alpha;Energy/keV;Counts",300,1700,2800);
	//能量刻度
	Double_t Cal_Parameter[3];
	TFile f_in1("C:/root/Projects/LaBr3/SGB_2inch_background_WithShield/Tree/EnergyCalibration.root");
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);
	f_in1.Close();
	//打开CCM的root文件，以获取部分电荷和总电荷
	TFile f_in("C:/root/Projects/LaBr3/SGB_2inch_background_WithShield/Tree/CCM.root");
	TTree *t_CCM = (TTree*) f_in.GetObjectChecked("t_CCM","TTree");
	Long64_t nentries = t_CCM->GetEntries();
	Float_t Q_total,Q_part;
	t_CCM->SetBranchAddress("Q_part",&Q_part);
	t_CCM->SetBranchAddress("Q_total",&Q_total);
	Long64_t nentries = t_CCM->GetEntries();
	Long64_t iwaveform;
	Int_t gamma_cnt=0,alpha_cnt=0;
	for(iwaveform=0;iwaveform < nentries;iwaveform++)
	{
		t_CCM->GetEntry(iwaveform);
		Float_t energy=Cal_Parameter[0]+Cal_Parameter[1]*Q_total+Cal_Parameter[2]*Q_total*Q_total;
		Float_t CCM_temp = Q_part/Q_total;//CCM比值作为甄别参数
		if (CCM_temp < CCM_thresh && energy>1600)
		{
			h_Spectrum_gamma->Fill(energy);
		}
		else if(energy>1600)
		{
			h_Spectrum_alpha->Fill(energy);
		}
	}
	cout<<nentries<<endl;
	gROOT->cd();
	gStyle->SetOptStat(kFALSE);
	//绘制甄别之后的gamma谱
	TCanvas *c_Gamma = new TCanvas("c_Sepration","CCM Distinction",1000,600);
//	c_Gamma->SetLogy();
	h_Spectrum_gamma->SetLineColor(kBlue);
	h_Spectrum_gamma->Draw("same");
	h_Spectrum_gamma->GetXaxis()->CenterTitle(kTRUE); 
	h_Spectrum_gamma->GetYaxis()->CenterTitle(kTRUE);
	//h_Spectrum_gamma->GetYaxis()->SetRangeUser(1,8000);
	//绘制甄别之后的alpha谱
	TCanvas *c_Alpha= new TCanvas("c_Actual","CCM Actual",1000,600);
//	c_Alpha->SetLogy();
	h_Spectrum_alpha->SetLineColor(kRed);
	h_Spectrum_alpha->Draw("same");
	//h_Spectrum_alpha->GetYaxis()->SetRangeUser(1,8000);
	h_Spectrum_alpha->GetXaxis()->CenterTitle(kTRUE); 
	h_Spectrum_alpha->GetYaxis()->CenterTitle(kTRUE);
	c_Gamma->SaveAs("C:/root/Projects/LaBr3/SGB_2inch_background_WithShield/Picture/Gamma.png");
	c_Alpha->SaveAs("C:/root/Projects/LaBr3/SGB_2inch_background_WithShield/Picture/Alpha.png");
	f_in.Close();

}
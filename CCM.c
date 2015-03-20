void CCM(char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield")
{
	gROOT->Reset();
	//全局变量写在前面，否则会因工作环境的切换导致图片画不出来！！！
	Float_t AlinedWave[400];
	Float_t Time[400];
	TH2F *h2D_CCMRate= new TH2F("h2D_CCMRate","TotalCharge VS CCM;TotalCharge;CCM",100,0,3000,30,NULL,NULL);
	//能量刻度
	Double_t Cal_Parameter[3];
	TFile f_in1(Form("%s/Tree/EnergyCalibration.root",filepath));
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);
	f_in1.Close();
	//选择部分电荷区间
	Float_t Time1=20,Time2=54.8;
	Float_t Deltatime = 0.4;
	TFile f_in(Form("%s/Tree/PreProcessedData.root",filepath));
	TTree *t1 = (TTree*) f_in.GetObjectChecked("t1","TTree");
	TBranch *b_time = t1->GetBranch("Time");
	b_time->SetAddress(Time);
	TBranch *b_ampl = t1->GetBranch("AlinedWave");
	b_ampl->SetAddress(AlinedWave);
	Float_t TotalCharge;
	t1->SetBranchAddress("TotalCharge",&TotalCharge);
	gStyle->SetOptStat(kFALSE);
	Long64_t iwaveform;
	Long64_t nentries = t1->GetEntries();
	Float_t Q_part,Q_total;
	TFile *f_CCMTree = new TFile(Form("%s/Tree/CCM.root",filepath),"recreate");
	TTree *t_CCM = new TTree("t_CCM","CCM data");
	t_CCM->Branch("Q_part",&Q_part,"Q_part/F");
	t_CCM->Branch("Q_total",&Q_total,"Q_total/F");
	for(iwaveform=0;iwaveform < nentries;iwaveform++)
	{
		t1->GetEntry(iwaveform);
		Float_t energy=Cal_Parameter[0]+Cal_Parameter[1]*TotalCharge+Cal_Parameter[2]*TotalCharge*TotalCharge;
		Float_t PartCharge = 0.0;
		for (Int_t i=(Int_t)(Time1/Deltatime);i<(Int_t)(Time2/Deltatime);i++)
		{
			PartCharge += AlinedWave[i]*0.4;
		}
		Float_t CCM_rate = PartCharge/TotalCharge;
		h2D_CCMRate->Fill(energy,CCM_rate);
		Q_part = PartCharge;
		Q_total = TotalCharge;
		t_CCM->Fill();
		if(iwaveform%2000==0) cout<<iwaveform<<"/"<<nentries<<endl;
	}
	t_CCM->Write();
 	gROOT->cd();
 	gStyle->SetOptStat(kFALSE);
	TCanvas *c3 = new TCanvas("c3","CCM RATE",1000,600);
	h2D_CCMRate->SetMarkerSize(0.8);
	h2D_CCMRate->GetXaxis()->CenterTitle(kTRUE); 
	h2D_CCMRate->GetYaxis()->CenterTitle(kTRUE); 
	h2D_CCMRate->Draw();
	c3->SaveAs(Form("%s/Picture/CCM Distinction.png",filepath));
	f_CCMTree->Close();
	f_in.Close();
}
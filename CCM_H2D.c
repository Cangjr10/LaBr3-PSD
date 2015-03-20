void CCM_H2D()
{
	gROOT->Reset();
	struct ROI_s 
	{
		Float_t LeftBoard;
		Float_t RightBoard;
		Float_t Center;
	};
	//��ȡ�����̶�����
	Double_t Cal_Parameter[3];
	TFile f_in1("C:/root/Projects/LaBr3/SGB_2inch/Tree/EnergyCalibration.root");
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);//��仰�ǳ���Ҫ��������ȡ��������������ӳ���ǲ����ģ�����ȷ��ʵ��

	const Int_t GammaPart1 = 100;//200-1000keV
	const Int_t GammaPart2=0;//1000-1200,1200-1300,1400-1500.
	const Int_t AlphaPart1=50;//1800-2700,
	const Int_t nTotalParts = GammaPart1/*+GammaPart2*/+AlphaPart1;
	ROI_s ROI[GammaPart1/*+GammaPart2*/+AlphaPart1];
	for (Int_t i=0;i<GammaPart1;i++)
	{
		ROI[i].LeftBoard = 300+i*(1800-300)/GammaPart1;
		ROI[i].RightBoard = 300+(i+1)*(1800-300)/GammaPart1;
		ROI[i].Center = (ROI[i].LeftBoard + ROI[i].RightBoard)/2;
	}
	for (Int_t i=GammaPart1;i<GammaPart1+AlphaPart1;i++)
	{
		ROI[i].LeftBoard = 1800+(i-(GammaPart1))*(2600-1800)/AlphaPart1 ;
		ROI[i].RightBoard = 1800+(i+1-(GammaPart1))*(2600-1800)/AlphaPart1 ;
		ROI[i].Center = (ROI[i].LeftBoard + ROI[i].RightBoard)/2;
	}
	//���仮�����
/******************************************************�Ե�ַ���ݽ�����ϣ��Եõ���˹��λ****************************************************************/
	TFile f_in("C:/root/Projects/LaBr3/SGB_2inch/Tree/CCM.root");
	TTree *t_CCM = (TTree*) f_in.GetObjectChecked("t_CCM","TTree");
	Long64_t nentries = t_CCM->GetEntries();
	Float_t Q_total,Q_part;
	t_CCM->SetBranchAddress("Q_part",&Q_part);
	t_CCM->SetBranchAddress("Q_total",&Q_total);
	Float_t CCM_tree;
	Float_t energy_tree;
	//�����treeֻ��Ϊ�˻�ͼ���㣬����CCM vs. Energy�����ӷֲ�
	TFile *f_out = new TFile("C:/root/Projects/LaBr3/SGB_2inch/Tree/CCMAnalysis.root","recreate");
	TTree *t_CCMAnalysis = new TTree("t_CCMAnalysis","tree of CCM Analysis");
	t_CCMAnalysis->Branch("CCM_tree",&CCM_tree,"CCM_tree/F");
	t_CCMAnalysis->Branch("energy_tree",&energy_tree,"energy_tree/F");
	Long64_t nentries = t_CCM->GetEntries();
	Long64_t iwaveform;
	for (Int_t i=0;i<nTotalParts;i++)
	{
		Float_t cnt=0;
		for(iwaveform=0;iwaveform < nentries;iwaveform++)
		{
			t_CCM->GetEntry(iwaveform);
			Float_t Energy0 = Cal_Parameter[0]+Cal_Parameter[1]*Q_total+Cal_Parameter[2]*Q_total*Q_total;
			if (Energy0 < ROI[i].RightBoard && Energy0 > ROI[i].LeftBoard)
			{
				if(cnt>200) continue;//��ÿ�������ֻ��ȡ���ּ�������ʹ�������ļ������ȣ��õ����ӵķֲ����ɣ��м�ƫ�죬����ƫ���ĸ�˹�ֲ���
				energy_tree = Energy0;
				CCM_tree = Q_part/Q_total;
				t_CCMAnalysis->Fill();
				cnt++;
			}
		}
		cout<<i<<"/"<<nTotalParts<<endl;
	}
	t_CCMAnalysis->Write();
	gROOT->cd();
	//���ƶ�ά���ӳ�����ľ��ȷֲ�
	TCanvas *c_CCM_dis = new TCanvas("c_CCM_dis","Distribution of CCM",960,0,1000,600);
	TH2F *h_2D = new TH2F("h_2D","CCM vs. Energy;Energy/keV;CCM",200,100,2600,100,0.68,0.72);
	h_2D->SetStats(kFALSE);
	h_2D->SetMaximum(16);
	h_2D->GetXaxis()->CenterTitle(kTRUE); 
	h_2D->GetYaxis()->CenterTitle(kTRUE);
	h_2D->GetYaxis()->SetTitleOffset(1.2);
	t_CCMAnalysis->Draw("CCM_tree:energy_tree>>h_2D","","COLZ");//CONT,COLZ
	c_CCM_dis->SaveAs("C:/root/Projects/LaBr3/SGB_2inch/Picture/CCM_Averaged.png");
	f_out->Close();
	f_in1.Close();
}
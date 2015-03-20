void CCM_GaussFit()
{
	gROOT->Reset();
	const Int_t Total_Bin_Counts = 50;
	TH1F *h_gamma = new TH1F("h_gamma","Charge Comparison Method;CCM_RATE;counts",Total_Bin_Counts,/*NULL,NULL*/0.685,0.725);//Time1=94.6,Time2=130;0.42,0.54
	TH1F *h_alpha = new TH1F("h_alpha","Charge Comparison Method;CCM_RATE;counts",Total_Bin_Counts,/*NULL,NULL*/0.685,0.725);//Time1=94.6,Time2=130;
	Bool_t gamma_flag;
	Float_t thresh = 0.7057;//甄别域以计算甄别率！
	//能量刻度
	Double_t Cal_Parameter[3];
	TFile f_in1("C:/root/Projects/LaBr3/SGB_2inch/Tree/EnergyCalibration.root");
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);
	f_in1.Close();
	TFile f_in("C:/root/Projects/LaBr3/SGB_2inch/Tree/CCM.root");
	TTree *t_CCM = (TTree*) f_in.GetObjectChecked("t_CCM","TTree");
	Long64_t nentries = t_CCM->GetEntries();
	Float_t Q_total,Q_part;
	t_CCM->SetBranchAddress("Q_part",&Q_part);
	t_CCM->SetBranchAddress("Q_total",&Q_total);
	Long64_t nentries = t_CCM->GetEntries();
	Long64_t iwaveform;
	Int_t gamma_cnt=0,alpha_cnt=0;
	Int_t gamma_right=0,alpha_right=0;//正确识别的alpha、gamma数目
	const Int_t wavenumber = 1999;
	for(iwaveform=0;iwaveform < nentries;iwaveform++)
	{
		t_CCM->GetEntry(iwaveform);
		Float_t energy=Cal_Parameter[0]+Cal_Parameter[1]*Q_total+Cal_Parameter[2]*Q_total*Q_total;
		if (energy<700 && energy>620)//energy<1500 && energy>1400
		{
			gamma_flag = true;
		}
		else if (energy<2600 && energy>1800)
		{
			gamma_flag = false;
		}
		else
		{
			continue;
		}
		if (gamma_flag)
		{
			if (gamma_cnt>wavenumber) continue;
			h_gamma->Fill(Q_part/Q_total);
			gamma_cnt++;
			if (Q_part/Q_total < thresh)
			{
				gamma_right++;
			}
		}
		else
		{
			if (alpha_cnt>wavenumber) continue;
			h_alpha->Fill(Q_part/Q_total);
			alpha_cnt++;
			if (Q_part/Q_total > thresh)
			{
				alpha_right++;
			}
		}
	}
	cout<<"gamma_eff = "<< gamma_right*10000/wavenumber<<endl;
	cout<<"alpha_eff = "<< alpha_right*10000/wavenumber<<endl;
 	gROOT->cd();
	gStyle->SetOptStat(kFALSE);
 	TCanvas *c1 = new TCanvas("CCM_Painter","CCM",960,0,1000,600);
	//寻找两张直方图的最大值，以设置合适的坐标范围
	Float_t max_axis1 = h_alpha->GetMaximum();
	Float_t max_axis2 = h_gamma->GetMaximum();
	Float_t max_axis = max_axis1>max_axis2 ? max_axis1 : max_axis2;
	h_gamma->SetLineColor(kBlack);
	h_gamma->SetLineStyle(7);
	h_gamma->SetLineWidth(2);
	h_gamma->GetYaxis()->SetRangeUser(0,1.1*max_axis);
	h_gamma->GetXaxis()->CenterTitle(kTRUE); 
	h_gamma->GetYaxis()->CenterTitle(kTRUE);
	h_gamma->GetXaxis()->SetTitleSize(0.045); 
	h_gamma->GetYaxis()->SetTitleSize(0.045); 
	h_gamma->Draw();
	h_alpha->SetLineColor(kBlack);
	h_alpha->SetLineStyle(1);
	h_alpha->SetLineWidth(2);
	h_alpha->GetXaxis()->CenterTitle(kTRUE); 
	h_alpha->GetYaxis()->CenterTitle(kTRUE);
	h_alpha->Draw("same");
	
	//绘制标识
	Int_t reference_bin = Int_t(h_gamma->GetMaximumBin()-Total_Bin_Counts*1/8);
	Float_t start_x = h_gamma->GetBinCenter(Int_t(reference_bin-Total_Bin_Counts*1/15));
	Float_t end_x = h_gamma->GetBinCenter(Int_t(reference_bin));
	Float_t start_y = 0.5*h_gamma->GetMaximum();
	Float_t end_y = start_y-h_gamma->GetMaximum()/10;
	TArrow ar1(start_x,start_y,end_x,end_y,0.02,"|>");
	ar1.SetLineColor(1);
	ar1.SetLineWidth(2);
	ar1.DrawClone();
	Float_t text_x = start_x - h_gamma->GetBinWidth(1)*Total_Bin_Counts*0.1;
	TLatex text(text_x,start_y,"Gamma Events");
	text.SetTextSize(0.04);
	text.DrawClone();
	Int_t reference_bin1 = Int_t(h_alpha->GetMaximumBin()+Total_Bin_Counts*1/8);
	Float_t start_x1 = h_alpha->GetBinCenter(Int_t(reference_bin1+Total_Bin_Counts*1/15));
	Float_t end_x1 = h_alpha->GetBinCenter(Int_t(reference_bin1));
	Float_t start_y1 = 0.7*h_gamma->GetMaximum();
	Float_t end_y1 = start_y1-h_gamma->GetMaximum()/10;
	Float_t text_x1 = start_x1 - h_alpha->GetBinWidth(1)*Total_Bin_Counts*0.08;
	TArrow ar2(start_x1,start_y1,end_x1,end_y1,0.02,"|>");
	ar2.SetLineColor(1);
	ar2.SetLineWidth(2);
	ar2.DrawClone();
	TLatex text2(text_x1,start_y1,"Alpha Events");
	text2.SetTextSize(0.04);
	text2.DrawClone();
	//甄别线;
	TArrow Arrow_thresh(thresh,max_axis*0.8,thresh,0,0.01,"|>");
	Arrow_thresh.SetLineColor(1);
	Arrow_thresh.SetLineWidth(1);
	Arrow_thresh.DrawClone();
	TLine* gLine1 = new TLine(thresh,max_axis*0.8,thresh-h_alpha->GetBinWidth(1),max_axis*0.8);
	gLine1->SetLineColor(1);
	gLine1->SetLineStyle(1);
	gLine1->SetLineWidth(1);
	gLine1->Draw();
	TLatex text1(thresh-h_alpha->GetBinWidth(1)*0.13*Total_Bin_Counts,max_axis*0.78,"Thresh");
	text1.SetTextSize(0.04);
	text1.DrawClone();
	c1->SaveAs("C:/root/Projects/LaBr3/SGB_2inch/Picture/CCM_Gauss.png");	
	//对数据进行高斯拟合
// 	TF1* g_gamma = new TF1("g_gamma","gaus"/*,NULL,NULL*/);
// 	TF1* g_alpha = new TF1("g_alpha","gaus"/*,NULL,NULL*/);
// 	h_gamma->Fit(g_gamma/*,"R"*/);
// 	h_alpha->Fit(g_alpha/*,"R"*/);
// 	Float_t FitParm_gammma[2];
// 	Float_t FitParm_alpha[2];
// 	//获取中心值
// 	FitParm_gammma[0] = g_gamma->GetParameter(1);
// 	FitParm_alpha[0]=g_alpha->GetParameter(1);
// 	//获取sigma
// 	FitParm_gammma[1] = g_gamma->GetParameter(2);
// 	FitParm_alpha[1]=g_alpha->GetParameter(2);
// 	Float_t FOM= (FitParm_alpha[0]-FitParm_gammma[0])/(2.355*(FitParm_alpha[1]+FitParm_gammma[1]));
// 	cout<<"FOM="<<TMath::Abs(FOM)<<endl;
	f_in.Close();
}
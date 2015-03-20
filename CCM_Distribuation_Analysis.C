#include <cmath>
Double_t fCCM(Double_t *x, Double_t *par)
{
	double xx=x[0];
	double p0=par[0];
	double p1=par[1];
	double p2=par[2];
	return p0+p1*xx+p2*xx*xx;
}
Double_t fSigma(Double_t *x, Double_t *par)
{
	double xx=x[0];
	double p0=par[0];
	double p1 = par[1];
	double p2=par[2];
	double rho=par[3];
	//double time_jitter = par[4];
	double CCM_Value = p0+p1*xx+p2*xx*xx;
	double sigma_rate = 0.314926/TMath::Sqrt(xx);//sigma Q/Q
	double statistical = TMath::Sqrt(2*(1-rho))*CCM_Value*sigma_rate;
	double noise_contribute = par[4]/xx;
	return  TMath::Sqrt(pow(statistical,2)+ pow(noise_contribute,2));
}

void CCM_Distribuation_Analysis(char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield",Float_t CCM_1 = 0.52,Float_t CCM_2=0.59)
{
	gROOT->Reset();
	struct ROI_s
	{
		Float_t LeftBoard;
		Float_t RightBoard;
		Float_t Center;
	};
	const Int_t GammaPart1 = 10;//300-1300keV
	const Int_t GammaPart2=1;//1400-1500.
	const Int_t AlphaPart1=7;//1800-2600,
	const Int_t nTotalParts = GammaPart1+GammaPart2+AlphaPart1;
	ROI_s ROI[GammaPart1+GammaPart2+AlphaPart1];
	for (Int_t i=0;i<GammaPart1;i++)
	{
		ROI[i].LeftBoard = 300+i*(1300-300)/GammaPart1;
		ROI[i].RightBoard = 300+(i+1)*(1300-300)/GammaPart1;
		ROI[i].Center = (ROI[i].LeftBoard + ROI[i].RightBoard)/2;
	}
	ROI[GammaPart1].LeftBoard = 1400;
	ROI[GammaPart1].RightBoard = 1500;
	ROI[GammaPart1].Center = 1450;
	// 	ROI[GammaPart1+1].LeftBoard = 1400;
	// 	ROI[GammaPart1+1].RightBoard = 1500;
	// 	ROI[GammaPart1+1].Center = 1450;
	// 	ROI[GammaPart1+1].LeftBoard = 1200;
	// 	ROI[GammaPart1+1].RightBoard = 1300;
	// 	ROI[GammaPart1+1].Center = 1250;
	// 	ROI[GammaPart1+2].LeftBoard = 1400;
	// 	ROI[GammaPart1+2].RightBoard = 1500;
	// 	ROI[GammaPart1+2].Center = 1450;
	for (Int_t i=GammaPart1+GammaPart2;i<GammaPart1+GammaPart2+AlphaPart1;i++)
	{
		ROI[i].LeftBoard = 1800+(i-(GammaPart1+GammaPart2))*(2500-1800)/AlphaPart1 ;
		ROI[i].RightBoard = 1800+(i+1-(GammaPart1+GammaPart2))*(2500-1800)/AlphaPart1 ;
		ROI[i].Center = (ROI[i].LeftBoard + ROI[i].RightBoard)/2;
	}
	//区间划分完毕
	Float_t mu[nTotalParts],sigma[nTotalParts],Energy[nTotalParts];
	Float_t Error_mu[nTotalParts],Error_sigma[nTotalParts];
	//能量刻度
	Double_t Cal_Parameter[3];
	TFile f_in1(Form("%s/Tree/EnergyCalibration.root",filepath));
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);
	TFile f_in(Form("%s/Tree/CCM.root",filepath));
	TTree *t_CCM = (TTree*) f_in.GetObjectChecked("t_CCM","TTree");
	Long64_t nentries = t_CCM->GetEntries();
	Float_t Q_total,Q_part;
	t_CCM->SetBranchAddress("Q_part",&Q_part);
	t_CCM->SetBranchAddress("Q_total",&Q_total);
	// 	Float_t CCM_tree;
	// 	Float_t energy_tree;
	// 	TFile *f_out = new TFile(Form("%s/Tree/CCMAnalysis.root","recreate");//Peak.root存储基线及幅度信息
	// 	TTree *t_CCMAnalysis = new TTree("t_CCMAnalysis","tree of CCM Analysis");
	// 	/*	t1->Branch("iwaveform",&iwaveform,"iwaveform/l");*/
	// 	t_CCMAnalysis->Branch("CCM_tree",&CCM_tree,"CCM_tree/F");
	// 	t_CCMAnalysis->Branch("energy_tree",&energy_tree,"energy_tree/F");
	Long64_t nentries = t_CCM->GetEntries();
	Long64_t iwaveform;
	for (Int_t i=0;i<nTotalParts;i++)
	{
// 		switch (i)
// 		{
// 		case 0:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.675,0.725);break;
// 		case 1:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.68,0.72);break;
// 		case 2:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.68,0.72);break;
// 		case 3:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.68,0.72);break;
// 		case 4:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.68,0.72);break;
// 		case 5:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.685,0.71);break;
// 		case 6:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.685,0.71);break;
// 		case 7:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.685,0.71);break;
// 		case 8:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.685,0.71);break;
// 		case 9:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.685,0.71);break;
// 		case 10:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.687,0.707);break;
// 			//alpha part
// 		case 11:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.72);break;
// 		case 12:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.72);break;
// 		case 13:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.72);break;
// 		case 14:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.72);break;
// 		case 15:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.715);break;
// 		case 16:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.715);break;
// 		default:
// 			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,0.695,0.715);
// 			break;
// 		}
		if (i < 11)
		{
			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,CCM_1,CCM_2);
		}
		else
		{
			TH1F *h_CCM = new TH1F("h_CCM","CCM Distinction",25,CCM_1,CCM_2);
		}
		
		Float_t cnt=0;
		for(iwaveform=0;iwaveform < nentries;iwaveform++)
		{
			t_CCM->GetEntry(iwaveform);
			Float_t Energy0 = Cal_Parameter[0]+Cal_Parameter[1]*Q_total+Cal_Parameter[2]*Q_total*Q_total;
			if (Energy0 < ROI[i].RightBoard && Energy0 > ROI[i].LeftBoard)
			{
				if(cnt>2499) continue;
				h_CCM->Fill(Q_part/Q_total);
				cnt++;
			}
		}
		cout<<cnt<<endl;
		TF1* g1 = new TF1("m1","gaus"/*,0.68,0.72*/);

		h_CCM->Fit(g1,"R");
		mu[i]=g1->GetParameter(1);
		Error_mu[i]=g1->GetParError(1);
		sigma[i]=g1->GetParameter(2);
		Error_sigma[i]=g1->GetParError(2);
		Energy[i] = ROI[i].Center;
		cout<<"Error_mu"<<Error_mu[i]<<endl;
		// 		if (i == 12)
		// 		{
		// 			gROOT->cd();
		// 			TCanvas *c_test = new TCanvas("CCM_test","CCM",1000,600);
		// 			h_CCM->Draw();
		// 			c_test->SaveAs(Form("%s/Picture/Gauss_test.png");
		// 			break;
		// 		}
	}
	//gROOT->cd();
	TF1 *f_mu_gamma = new TF1("f_mu_gamma",fCCM,300,1500,3);
	TF1 *f_mu_alpha = new TF1("f_mu_alpha",fCCM,1800,2500,3);
	//TCanvas *c_mu = new TCanvas("CCM mu Distribution","CCM_mu",960,0,960,600);
	TGraphErrors *g_mu = new TGraphErrors(nTotalParts,Energy,mu,NULL,Error_mu);
	f_mu_gamma->SetParameter(0,0.65);
	f_mu_gamma->SetParameter(1,0);
	f_mu_gamma->SetParameter(2,0);
	g_mu->Fit("f_mu_gamma","RN");
	double Parameter_mu_Gamma[3];
	f_mu_gamma->GetParameters(&Parameter_mu_Gamma[0]);//获取拟合参数
	f_mu_alpha->SetParameter(0,CCM_1);
	f_mu_alpha->SetParameter(1,0);
	f_mu_alpha->SetParameter(2,0);
	f_mu_alpha->SetParLimits(2,-1E-8,0);
	g_mu->Fit("f_mu_alpha","RN");
	
	double Parameter_mu_Alpha[3];
	f_mu_alpha->GetParameters(&Parameter_mu_Alpha[0]);//获取拟合参数
	//gROOT->Reset();
	gROOT->cd();
	TCanvas *c_mu = new TCanvas("CCM mu Distribution","CCM_mu",0,0,1440,1080);
	//画图部分
	g_mu->SetTitle("Mean Value vs. Energy for CCM;Energy/keV;Mean Value");
	g_mu->SetLineColor(kBlack);
	g_mu->SetMarkerColor(kBlack);
	g_mu->SetMarkerSize(1);
	//g_mu->SetMarkerStyle(21);
	g_mu->Draw("APE*");
	g_mu->GetXaxis()->CenterTitle(kTRUE); 
	g_mu->GetYaxis()->CenterTitle(kTRUE);
	g_mu->GetYaxis()->SetRangeUser(CCM_1,CCM_2);
	g_mu->GetYaxis()->SetTitleOffset(1.2);
	//补画gamma部分的拟合曲线
	TF1 *f_mu_gamma1 = new TF1("f_mu_gamma1",fCCM,300,2500,3);
	for (Int_t i=0;i<3;i++)
	{
		f_mu_gamma1->SetParameter(i,Parameter_mu_Gamma[i]);
	}
	f_mu_gamma1->SetLineColor(kBlack);
	f_mu_gamma1->SetLineStyle(1);
	f_mu_gamma1->SetLineWidth(1);
	f_mu_gamma1->Draw("same");
// 	TF1 *f_mu_gamma2 = new TF1("f_mu_gamma1",fCCM,1500,2500,3);
// 	for (Int_t i=0;i<3;i++)
// 	{
// 		f_mu_gamma2->SetParameter(i,Parameter_mu_Gamma[i]);
// 	}
// 	f_mu_gamma2->SetLineColor(kBlack);
// 	f_mu_gamma2->SetLineWidth(2);
// 	f_mu_gamma2->SetLineStyle(1);
// 	f_mu_gamma2->Draw("same");
	TF1 *f_mu_alpha1 = new TF1("f_mu_alpha1",fCCM,1800,2500,3);
	for (Int_t i=0;i<3;i++)
	{
		f_mu_alpha1->SetParameter(i,Parameter_mu_Alpha[i]);
	}
	f_mu_alpha1->SetLineColor(kBlack);
	f_mu_alpha1->SetLineWidth(2);
	f_mu_alpha1->SetLineStyle(7);
	f_mu_alpha1->Draw("same");
	TLegend *legend=new TLegend(0.75,0.78,0.90,0.88);//
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->AddEntry(f_mu_alpha1,"alpha","l");
	legend->AddEntry(f_mu_gamma1,"gamma","l");
	legend->Draw();
	c_mu->SaveAs(Form("%s/Picture/Mu Distribution.pdf",filepath));

// 	/****************************************************sigma部分***********************************************************************/
	TCanvas *c_sigma = new TCanvas("CCM sigma Distribution","CCM_sigma",0,0,1440,1080);
	double Parameter_sigma_Gamma[5];
	double Parameter_sigma_Alpha[5];
	//alpha,gamma的拟合函数
	TF1 *f_sigma_gamma = new TF1("f_sigma_gamma",fSigma,300,1500,5);
	TF1 *f_sigma_alpha = new TF1("f_sigma_alpha",fSigma,1800,2500,5);
	/*	TGraph *g_sigma = new TGraph(nTotalParts,Energy,sigma);*/
	TGraphErrors *g_sigma = new TGraphErrors(nTotalParts,Energy,sigma,NULL,Error_sigma);
	//对gamma段的sigma分布进行拟合，带入中心值
	f_sigma_gamma->FixParameter(0,Parameter_mu_Gamma[0]);
	f_sigma_gamma->FixParameter(1,Parameter_mu_Gamma[1]);
	f_sigma_gamma->FixParameter(2,Parameter_mu_Gamma[2]);
	f_sigma_gamma->SetParameter(3,0.9);
	//f_sigma_gamma->FixParameter(4,0);
	g_sigma->Fit("f_sigma_gamma","RN","ep");
	f_sigma_gamma->GetParameters(&Parameter_sigma_Gamma[0]);//获取拟合参数
	//alpha拟合，带入中心值
	f_sigma_alpha->FixParameter(0,Parameter_mu_Alpha[0]);
	f_sigma_alpha->FixParameter(1,Parameter_mu_Alpha[1]);
	f_sigma_alpha->FixParameter(2,Parameter_mu_Alpha[2]);
	f_sigma_alpha->SetParameter(3,0.9);
	//f_sigma_alpha->FixParameter(4,0);
	g_sigma->Fit("f_sigma_alpha","RN","ep");
	f_sigma_alpha->GetParameters(&Parameter_sigma_Alpha[0]);//获取拟合参数

	//绘图
	g_sigma->Draw("APE*");
	g_sigma->SetTitle("Sigma vs. Energy for CCM;Energy/keV;Sigma");
	g_sigma->SetLineColor(kBlack);
	g_sigma->SetMarkerColor(kBlack);
	g_sigma->SetMarkerSize(1);
	g_sigma->GetXaxis()->CenterTitle(kTRUE); 
	g_sigma->GetYaxis()->CenterTitle(kTRUE);
	//g_sigma->GetYaxis()->SetRangeUser(0.0015,0.006);
	g_sigma->GetYaxis()->SetTitleOffset(1.6);
	//补画gamma部分的拟合曲线
	TF1 *f_sigma_gamma_temp = new TF1("f_sigma_gamma_temp",fSigma,300,2500,5);
	/*	Parameter_sigma_Gamma[2]=0.1213;*/
	for (Int_t i=0;i<5;i++)
	{
		f_sigma_gamma_temp->SetParameter(i,Parameter_sigma_Gamma[i]);
	}
	f_sigma_gamma_temp->SetLineColor(kBlack);
	f_sigma_gamma_temp->SetLineWidth(2);
	f_sigma_gamma_temp->SetLineStyle(1);
	f_sigma_gamma_temp->Draw("same");
	TF1 *f_sigma_alpha_temp = new TF1("f_sigma_alpha_temp",fSigma,1800,2500,5);
	/*	Parameter_sigma_Alpha[2]=0.1213;*/
	for (Int_t i=0;i<5;i++)
	{
		f_sigma_alpha_temp->SetParameter(i,Parameter_sigma_Alpha[i]);
	}
	f_sigma_alpha_temp->SetLineColor(kBlack);
	f_sigma_alpha_temp->SetLineWidth(2);
	f_sigma_alpha_temp->SetLineStyle(7);
	f_sigma_alpha_temp->Draw("same");
	TLegend *legend=new TLegend(0.75,0.78,0.90,0.88);//
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->AddEntry(f_sigma_alpha_temp,"alpha","l");
	legend->AddEntry(f_sigma_gamma_temp,"gamma","l");
	legend->Draw();
 	c_sigma->SaveAs(Form("%s/Picture/Sigma Distribution.pdf",filepath));
	//绘制理想分布曲线图
	TCanvas *c_theory = new TCanvas("CCM theory Distribution","CCM Theory Distribution",0,0,1440,1080);
	TGraphErrors *g_ccm = new TGraphErrors(nTotalParts,Energy,mu,NULL,sigma);
	g_ccm->Draw("A*");
	g_ccm->GetYaxis()->SetRangeUser(CCM_1,CCM_2);
	const Int_t slice = 50;
	Float_t Energy_theory[slice];
	Float_t mu_thoery[slice];
	Float_t sigma_theory[slice];
	for (Int_t i=0;i<slice;i++)
	{
		Energy_theory[i] = 300+(i+1)*(2500-300)/slice;
		mu_thoery[i] = Parameter_mu_Gamma[0]+Parameter_mu_Gamma[1]*Energy_theory[i]+Parameter_mu_Gamma[2]*Energy_theory[i]*Energy_theory[i];
		sigma_theory[i] = TMath::Sqrt(pow(mu_thoery[i]*Parameter_sigma_Gamma[3]/TMath::Sqrt(Energy_theory[i]),2) + pow(Parameter_sigma_Gamma[4],2));
	}
	TGraphErrors *g_theory_gamma = new TGraphErrors(slice,Energy_theory,mu_thoery,NULL,sigma_theory);
	g_theory_gamma->SetFillColor(kBlue-4);
	g_theory_gamma->Draw("sameE3L");
	g_theory_gamma->SetTitle("Inferential Distribution of CCM vs.Energy;Energy/keV;CCM");
	//alpha 段理论分布

	for (Int_t i=0;i<slice;i++)
	{
		Energy_theory[i] = 1800+(i+1)*(2500-1800)/slice;
		mu_thoery[i] = Parameter_mu_Alpha[0]+Parameter_mu_Alpha[1]*Energy_theory[i]+Parameter_mu_Alpha[2]*Energy_theory[i]*Energy_theory[i];
		sigma_theory[i] = TMath::Sqrt(pow(mu_thoery[i]*Parameter_sigma_Alpha[3]/TMath::Sqrt(Energy_theory[i]),2) + pow(Parameter_sigma_Alpha[4],2));
	}
	TGraphErrors *g_theory_alpha = new TGraphErrors(slice,Energy_theory,mu_thoery,NULL,sigma_theory);
	g_theory_alpha->SetFillColor(kRed-4);
	//g_theory_alpha->SetLineColor(kWhite);
	//g_theory_alpha->SetFillStyle(3004);
	g_theory_alpha->Draw("sameE3L");
	g_ccm->SetTitle("Inferential Distribution of CCM vs.Energy;Energy/keV;CCM");
	g_ccm->SetLineColor(kYellow-7);
	g_ccm->SetMarkerColor(kYellow-7);
	g_ccm->GetXaxis()->CenterTitle(kTRUE); 
	g_ccm->GetYaxis()->CenterTitle(kTRUE);
	g_ccm->GetYaxis()->SetTitleOffset(1.2);
	g_ccm->SetLineStyle(7);
	g_ccm->Draw("same*");
	f_mu_gamma1->SetLineColor(kWhite);
	f_mu_gamma1->Draw("sameA");
	f_mu_alpha1->SetLineColor(kWhite);
	f_mu_alpha1->Draw("sameA");
	c_theory->SaveAs(Form("%s/Picture/Theory Distribution.pdf",filepath));
	//	f_out->Close();
	f_in.Close();
	f_in1.Close();
}
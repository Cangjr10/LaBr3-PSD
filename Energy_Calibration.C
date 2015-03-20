Double_t fbkg(Double_t *x, Double_t *par)
{
	double xx=x[0];
	double p0=par[0];
	double p1=par[1];
	double p2=par[2];
	return p0+p1*xx+p2*xx*xx;
}

Double_t fsig(Double_t *x, Double_t *par)
{
	double xx=x[0];
	double mu=par[0];
	double sigma=par[1]/2.355;
	//double binwidth=9;
	double norm=par[2];

	return norm*TMath::Exp(-(xx-mu)*(xx-mu)/2/sigma/sigma)/TMath::Sqrt(2*TMath::Pi())/sigma;
}

Double_t fitf(Double_t *x, Double_t *par)
{
	return fbkg(x,par)+fsig(x,&par[3]);
}
void Energy_Calibration(Int_t calibration_arry_OK = 1,char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield")
{
	gROOT->Reset();
/*****************************************************************************如果能谱中不含高斯峰，给定刻度数组即可*****************************************/
	if (calibration_arry_OK == 1)//给定刻度的数组，只需进行二阶多项式拟合即可
	{
		const int ROI_Number = 3;//能量峰数目
		Float_t ROI_Energy_Centroid[ROI_Number] = {809,1765,2614.7};//能量
		Float_t ROI_Charge_Centroid[ROI_Number] = {2562,5555,8166};//能谱中观察出的电荷量的中心峰位，环境测量
		Double_t Cal_Parameter[3];//刻度曲线系数
		//读取preprocessdata.root绘制刻度之后的能谱图
		TH1F* h_Spectrum_energy = new TH1F("h_Spectrum_energy","SGB 2inch LaBr3;Energy/keV;Counts",1500,200,2800);
		TFile f_in(Form("%s/Tree/PreProcessedData.root",filepath));
		TTree *t1 = (TTree*) f_in.GetObjectChecked("t1","TTree");
		Float_t TotalCharge;
		t1->SetBranchAddress("TotalCharge",&TotalCharge);
		//保存能量刻度曲线
		TFile *f_out = new TFile(Form("%s/Tree/EnergyCalibration.root",filepath),"recreate");
		TTree *t_energy_cali = new TTree("t_energy_cali","energy_cali");
		t_energy_cali->Branch("ROI_Number",&ROI_Number,"ROI_Number/I");
		t_energy_cali->Branch("ROI_Energy_Centroid",ROI_Energy_Centroid,"ROI_Energy_Centroid[3]/F");
		t_energy_cali->Branch("ROI_Charge_Centroid",ROI_Charge_Centroid,"ROI_Charge_Centroid[3]/F");
		t_energy_cali->Branch("Cal_Parameter",Cal_Parameter,"Cal_Parameter[3]/D");
		//二阶多项式能量刻度
		TF1 f_EnergyCalib("POLY2 law","[0]+x*[1]+x*x*[2]",200,9000);
		TGraph *energy_cal_curve = new TGraph(ROI_Number,ROI_Charge_Centroid,ROI_Energy_Centroid);
		energy_cal_curve->Fit(&f_EnergyCalib);
		f_EnergyCalib.GetParameters(&Cal_Parameter[0]);
		t_energy_cali->Fill();
		t_energy_cali->Write();
		Long64_t iwaveform;
		Long64_t nentries = t1->GetEntries();
		for(iwaveform=0;iwaveform < nentries;iwaveform++)
		{
			t1->GetEntry(iwaveform);
			Float_t energy = Cal_Parameter[0]+Cal_Parameter[1]*TotalCharge+Cal_Parameter[2]*TotalCharge*TotalCharge;
			h_Spectrum_energy->Fill(energy);
			if(iwaveform%1000==0) cout << iwaveform << "/" << nentries <<endl;	
		}
		h_Spectrum_energy->Smooth(2);
		gROOT->cd();
		gStyle->SetOptStat(kFALSE);
		TCanvas* c_Spectrum_energy = new TCanvas("c_Spectrum_energy","Spectrum_energy Histo",0,0,1440,1080);
		h_Spectrum_energy->SetTitle("Intrinsic and Natural Background;Energy/keV;Counts");
		//c_Spectrum_energy->SetLogy();
		h_Spectrum_energy->GetXaxis()->CenterTitle(kTRUE); 
		h_Spectrum_energy->GetYaxis()->CenterTitle(kTRUE);
		h_Spectrum_energy->GetYaxis()->SetTitleOffset(1);
		h_Spectrum_energy->GetXaxis()->SetTitleSize(0.045); 
		h_Spectrum_energy->GetYaxis()->SetTitleSize(0.045); 
		h_Spectrum_energy->Draw();
		c_Spectrum_energy->SaveAs(Form("%s/Picture/NaturalBackground.pdf",filepath));
		f_in.Close();
		f_out->Close();
	}
/***********************************************************如果能谱中不含高斯峰，给定刻度数组即可****************************************************************/
	else//否则从电荷谱中，只需给出峰的大致中心位置即可，程序会根据拟合出来的中心值进行能量刻度
	{
		const int ROI_Number = 3;//能量峰数目
		Float_t ROI_Energy_Centroid[ROI_Number] = {512,661.6,1274.5};//
		Float_t ROI_Charge_Centroid_Ref[ROI_Number] = {1623.4,2094.9,4017.6};//能谱中观察出的电荷量的中心峰位
		Float_t ROI_Charge_Centroid[ROI_Number];//拟合得到的峰位中心
		Float_t Resolution_Ref = 0.03;//参考能量分辨率，以给出大致的拟合范围
		Float_t Factor_FWHM = 2;//区间选取宽度，左右各2.5FWHM
		TH1F* h_Spectrum = new TH1F("h_Spectrum","SGB 2inch LaBr3 Na-22 Cs-137;Charge;Counts",1000,0,9000);
		TFile f_in(Form("%s/Tree/PreProcessedData.root",filepath));
		TTree *t1 = (TTree*) f_in.GetObjectChecked("t1","TTree");
		Float_t TotalCharge;
		t1->SetBranchAddress("TotalCharge",&TotalCharge);
		Double_t Cal_Parameter[3];//刻度曲线系数
		TFile *f_out = new TFile(Form("%s/Tree/EnergyCalibration.root",filepath),"recreate");
		TTree *t_energy_cali = new TTree("t_energy_cali","energy_cali");
		t_energy_cali->Branch("ROI_Number",&ROI_Number,"ROI_Number/I");
		t_energy_cali->Branch("ROI_Energy_Centroid",ROI_Energy_Centroid,"ROI_Energy_Centroid[3]/F");
		t_energy_cali->Branch("ROI_Charge_Centroid",ROI_Charge_Centroid,"ROI_Charge_Centroid[3]/F");
		t_energy_cali->Branch("Cal_Parameter",Cal_Parameter,"Cal_Parameter[3]/D");
		Long64_t iwaveform;
		Long64_t nentries = t1->GetEntries();
		for(iwaveform=0;iwaveform < nentries;iwaveform++)
		{
			t1->GetEntry(iwaveform);
			h_Spectrum->Fill(TotalCharge);
		}
		//h_Spectrum->GetYaxis()->SetRangeUser(0,2300);
/******************************************************对道址数据进行拟合，以得到高斯峰位****************************************************************/
		TCanvas* c_Spectrum_energy = new TCanvas("c_Spectrum_energy","Spectrum_energy Histo",800,0,1000,600);
		Float_t ROI_Charge_left[ROI_Number];
		Float_t ROI_Charge_right[ROI_Number];
		Float_t ROI_Charge_FWHM[ROI_Number];
		for(Int_t i=0;i<ROI_Number;i++)
		{
			ROI_Charge_FWHM[i] = ROI_Charge_Centroid_Ref[i]*Resolution_Ref;
			ROI_Charge_left[i] = ROI_Charge_Centroid_Ref[i] - ROI_Charge_FWHM[i]*Factor_FWHM;
			ROI_Charge_right[i] = ROI_Charge_Centroid_Ref[i] + ROI_Charge_FWHM[i]*Factor_FWHM;
		}
		char* f_name[3]={"f1","f2","f3"};
		char* bkg_name[3]={"myfbkg1","myfbkg2","myfbkg3"};
		char* sig_name[3]={"myfsig1","myfsig2","myfsig3"};
		TF1 * myfbkg[ROI_Number];
		TF1 * myfsig[ROI_Number];
		TF1 * f_total[ROI_Number];
		for(Int_t i=0;i<ROI_Number;i++)
		{
			TF1 * f_total_temp = new TF1("channel",fitf,ROI_Charge_left[i],ROI_Charge_right[i],6);
			f_total_temp->SetParameter(3,ROI_Charge_Centroid_Ref[i]);
			f_total_temp->SetParameter(4,ROI_Charge_FWHM[i]);
			h_Spectrum->Fit("channel","R"/*,"ep"*/);
			double mypar[6];
			f_total_temp->GetParameters(&mypar[0]);
			ROI_Charge_Centroid[i] = mypar[3];
		}
/*************************************************************根据拟合得到高斯峰位进行能量刻度****************************************************************/

		const int n_points=ROI_Number;
		TF1 f_EnergyCalib("POLY2 law","[0]+x*[1]+x*x*[2]",200,9000);
		TGraph *energy_cal_curve = new TGraph(n_points,ROI_Charge_Centroid,ROI_Energy_Centroid);
		energy_cal_curve->Fit(&f_EnergyCalib);
		f_EnergyCalib.GetParameters(&Cal_Parameter[0]);
		t_energy_cali->Fill();
		t_energy_cali->Write();
		gROOT->cd();
		gStyle->SetOptStat(kFALSE);
		TH1F* h_Spectrum_energy = new TH1F("h_Spectrum_energy","SGB 2inch LaBr3 Co-60;Energy/keV;Counts",1000,0,3000);
		for(iwaveform=0;iwaveform < nentries;iwaveform++)
		{
			t1->GetEntry(iwaveform);
			Float_t energy = Cal_Parameter[0]+Cal_Parameter[1]*TotalCharge+Cal_Parameter[2]*TotalCharge*TotalCharge;
			h_Spectrum_energy->Fill(energy);
		}
		//能量域拟合获取特征参数，以计算能量分辨率
		Float_t ROI_Energy_left[ROI_Number];
		Float_t ROI_Energy_right[ROI_Number];
		Float_t ROI_Energy_FWHM[ROI_Number];
		for(Int_t i=0;i<ROI_Number;i++)
		{
			ROI_Energy_FWHM[i] = ROI_Energy_Centroid[i]*Resolution_Ref;
			ROI_Energy_left[i] = ROI_Energy_Centroid[i] - ROI_Energy_FWHM[i]*Factor_FWHM;
			ROI_Energy_right[i] = ROI_Energy_Centroid[i] + ROI_Energy_FWHM[i]*Factor_FWHM;
		}
		for(Int_t i=0;i<ROI_Number;i++)
		{
			f_total[i] = new TF1(f_name[i],fitf,ROI_Energy_left[i],ROI_Energy_right[i],6);
			f_total[i]->SetParameter(3,ROI_Energy_Centroid[i]);
			f_total[i]->SetParameter(4,ROI_Energy_FWHM[i]);
			h_Spectrum_energy->Fit(f_name[i],"R"/*,"ep"*/);
		}
		//c_Spectrum_energy->cd();
		h_Spectrum_energy->SetLineColor(6);
		h_Spectrum_energy->SetTitle("SGB 2inch LaBr3;Energy/keV;Counts");
		h_Spectrum_energy->GetXaxis()->CenterTitle(kTRUE); 
		h_Spectrum_energy->GetYaxis()->CenterTitle(kTRUE);
		h_Spectrum_energy->GetYaxis()->SetTitleOffset(1.2);
		h_Spectrum_energy->Draw("ep");
		for (Int_t i=0;i<ROI_Number;i++)
		{
			//获取拟合参数绘制本底
			double mypar[6];
			f_total[i]->GetParameters(&mypar[0]);
			cout<<"resolution="<<mypar[4]/mypar[3]<<endl;
			char res[100];
			sprintf(res,"resolution=%.1f%% @%.1f keV  ",100*mypar[4]/mypar[3],mypar[3]);  // float 到 char
			TLatex text(mypar[3]+100,h_Spectrum_energy->GetBinContent((Int_t)(mypar[3]*1000/3000)),res);//将电荷量转换至道数
			text.DrawClone();
			//绘制本底与信号曲线
			myfbkg[i] = new TF1(bkg_name[i],fbkg,ROI_Energy_left[i],ROI_Energy_right[i],3);
			myfsig[i] = new TF1(sig_name[i],fsig,ROI_Energy_left[i],ROI_Energy_right[i],3);
			for (int j=0;j<6;j++) 
			{
				f_total[i]->SetParameter(j,mypar[j]);
				if (j<3) myfbkg[i]->SetParameter(j,mypar[j]);
				else myfsig[i]->SetParameter(j-3,mypar[j]);
			}
			f_total[i]->Draw("same");
			f_total[i]->SetLineColor(2);
			myfbkg[i]->Draw("same");
			myfbkg[i]->SetLineColor(4);
		}
		c_Spectrum_energy->SaveAs(Form("%s/Picture/EnergySpectrum.png",filepath));
		f_out->Close();
		f_in.Close();
	}
	
}

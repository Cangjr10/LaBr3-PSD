void PreProcess_CFD(char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield") //电荷量归一化并将波形按最大值对齐
{
	gROOT->Reset();  
	struct data_s
	{
		Float_t time[1252];
		Float_t ampl[1252];
	};
	data_s data_in;
	TH1F *h_Spectrum = new TH1F("h_Spectrum","SGB 2inch LaBr3",1000,0,9000);
	cout<<filepath<<endl;
	TFile f_in(Form("%s/Tree/data.root",filepath));
	TTree *t0 = (TTree*) f_in.GetObjectChecked("t0","TTree");
	TBranch *b_data = t0->GetBranch("data");
	b_data->SetAddress(&data_in);
	Long64_t nentries = t0->GetEntries();
	//存储归一化之后的波形数据
	Long64_t iwaveform;
	Float_t CFD_value = 0.1;
	const Int_t WaveLength = 400;
	Float_t AlinedWave[400];
	Float_t AveragedWave[400];
	Float_t Time[400];
	Float_t DeltaTime = 0.4;
	Float_t TotalCharge = 1;
	Float_t AlineTime = 20;//ns
	Int_t start_loc;
	TFile *f_out = new TFile(Form("%s/Tree/PreProcessedData.root",filepath),"recreate");//Peak.root存储基线及幅度信息
	TTree *t1 = new TTree("t1","tree of smoothed preprocessed data");
/*	t1->Branch("iwaveform",&iwaveform,"iwaveform/l");*/
	t1->Branch("Time",Time,"Time[400]/F");
	t1->Branch("AlinedWave",AlinedWave,"AlinedWave[400]/F");
	t1->Branch("TotalCharge",&TotalCharge,"TotalCharge/F");
	t1->Branch("AveragedWave",AveragedWave,"AveragedWave[400]/F");
	//电荷量能谱图
	
	for(iwaveform=0;iwaveform < nentries;iwaveform++) 
	{
		t0->GetEntry(iwaveform);
		Float_t pedestal = TMath::Mean(200,data_in.ampl);//取前200个数，平均作为基线,因为触发不一样
		Float_t peakampl;
		peakampl = pedestal - (TMath::MinElement(1252,data_in.ampl));//计算幅度
		if (peakampl > 200 || peakampl<20) continue;
		//上升沿恒比对齐
		//寻找起始位置，20%幅度位置
		//记20%幅度处时间为40ns
		Int_t peak_loc = TMath::LocMin(1252,data_in.ampl);
// 		for (Int_t i = peak_loc; i>0; i--)
// 		{
// 			//寻找20% 幅度位置
// 			if((pedestal-data_in.ampl[i]) < CFD_value*peakampl)
// 			{
// 				start_loc = i;
// 				break;
// 			}
// 		}
		for (Int_t i = 0; i<peak_loc; i++)
		{
			//寻找20% 幅度位置
			if((pedestal-data_in.ampl[i]) > CFD_value*peakampl)
			{
				start_loc = i;
				break;
			}
		}

		Int_t Aline_loc = Int_t(AlineTime/DeltaTime);
		if (start_loc > Aline_loc && start_loc<1252-WaveLength+2)//防止数组越界
		{
			TotalCharge = 0.0;
			for (Int_t i=0;i<WaveLength;i++)//存储归一化幅度信息，变为正幅度
			{
				//以20%幅度处位置为100，即40ns
				AlinedWave[i] = pedestal-data_in.ampl[start_loc+i-Aline_loc];
				Time[i] = i*DeltaTime;
				TotalCharge += AlinedWave[i]*DeltaTime;
				AveragedWave[i] = AlinedWave[i]/peakampl;
			}
		}
 		h_Spectrum->Fill(TotalCharge);
		t1->Fill();
		if(iwaveform%1000==0) cout << iwaveform << "/" << nentries <</* " TotalCharge="<<TotalCharge<<*/" start_loc="<<start_loc<<" peak_loc="<<peak_loc<<endl;	
 	}
	t1->Write();
	cout << "Tree file generated." <<endl;

 	gROOT->cd();
 	TCanvas *c1 = new TCanvas("c1","A canvas",1000,600);
  	gStyle->SetOptStat(kFALSE);
	h_Spectrum->GetXaxis()->CenterTitle(kTRUE); 
	h_Spectrum->GetYaxis()->CenterTitle(kTRUE);
	h_Spectrum->GetYaxis()->SetTitleSize(0.045);
  	h_Spectrum->Draw();
	//波形累积
	/*gROOT->cd();*/
	TCanvas *c2 = new TCanvas("c2","accumulation of waveforms",960,0,960,600);
	TH2F *h_2D = new TH2F("h_2D","Accumulated Specified Normalized Waveforms;Time;Amplitude/mV",400,0,160,400,-10,200);
	h_2D->SetStats(kFALSE);
	h_2D->GetXaxis()->SetTitle("Time/ns");
	h_2D->GetYaxis()->SetTitle("Amplitude/mV");
	h_2D->GetXaxis()->CenterTitle(kTRUE); 
	h_2D->GetYaxis()->CenterTitle(kTRUE);
	h_2D->SetMaximum(200);
	t1->Draw("AlinedWave:Time>>h_2D","","COLZ");
	c1->SaveAs(Form("%s/Picture/Preprocess_CFD_Spectrum.png"),filepath);
	c2->SaveAs(Form("%s/Picture/Preprocess_CFD_Acculumation.png"),filepath);
	f_out->Close();
	f_in.Close();
}
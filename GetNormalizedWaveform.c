void GetNormalizedWaveform(char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield")
{
	gROOT->Reset();
	Float_t Time[400];
	Float_t AlinedWave[400];
	Float_t NormalizedTime[400];
	Float_t NormalizedWave[400];//��ʱ����
	Float_t TotalWave_alpha[400];
	Float_t TotalWave_gamma[400];//�ܵ�ɻ��ֲ���
	Float_t Average_alpha[400];//ƽ������
	Float_t Average_gamma[400];
	Float_t diff[400];//(alpha-gamma)/gamma;
	Int_t alpha_cnt=0;
	Int_t gamma_cnt=0;
	//��ȡ�����̶�����
	Double_t Cal_Parameter[3];
	TFile f_in1(Form("%s/Tree/EnergyCalibration.root",filepath));
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);//��仰�ǳ���Ҫ��������ȡ��������������ӳ���ǲ����ģ�����ȷ��ʵ��
	f_in1.Close();
// 	TH2F *h_2D = new TH2F("h_2D","Accumulated Specified Normalized Waveforms;Time;Amplitude/mV",400,0,240,400,0,150);
	memset(TotalWave_gamma,0,400);
	memset(TotalWave_alpha,0,400);
/******************************************************�����Ͳ������ı���ʽ���*****************************************/
	ofstream f_out_alpha(Form("%s/Text/wave_alpha.txt",filepath),ios_base::out);
	ofstream f_out_gamma(Form("%s/Text/wave_gamma.txt",filepath),ios_base::out);
	f_out_alpha <<setiosflags(ios_base::left)<<setw(10)<< "alpha event"<<endl;
	f_out_alpha <<setiosflags(ios_base::left)<<setw(10)<< "Time"<<resetiosflags(ios_base::left) <<setw(10)<<"Amp"<<endl;
	f_out_gamma << "gamma event"<<endl;
	f_out_gamma <<setiosflags(ios_base::left)<<setw(10)<< "Time"<<resetiosflags(ios_base::left) <<setw(10)<<"Amp"<<endl;
/******************************************************��ȡ�����ļ�*****************************************/
	TFile f_in(Form("%s/Tree/PreProcessedData.root",filepath));
	TTree *t1 = (TTree*) f_in.GetObjectChecked("t1","TTree");
	TBranch *b_time = t1->GetBranch("Time");
	b_time->SetAddress(Time);
	TBranch *b_ampl = t1->GetBranch("AlinedWave");
	b_ampl->SetAddress(AlinedWave);
	Float_t TotalCharge;
	t1->SetBranchAddress("TotalCharge",&TotalCharge);
	Int_t gamma_flag;
/******************************************************�����Ͳ��α���Ϊroot�ļ�������100��5��s-kƽ��*****************************************/
	TFile f_out(Form("%s/Tree/NormalizedData.root",filepath),"recreate");//Peak.root�洢���߼�������Ϣ
	TTree *t2 = new TTree("t2","tree of Normalized data");
	t2->Branch("NormalizedTime",NormalizedTime,"NormalizedTime[400]/F");
	t2->Branch("Average_alpha",Average_alpha,"Average_alpha[400]/F");
	t2->Branch("Average_gamma",Average_gamma,"Average_gamma[400]/F");
/*	gStyle->SetOptStat(kFALSE);*/
	Long64_t iwaveform;
	const Int_t wavenumber = 1999;//ȡ��ƽ�����ε���Ŀ
	Long64_t nentries = t1->GetEntries();
	for(iwaveform=0;iwaveform < nentries;iwaveform++)
	{
		if (alpha_cnt > wavenumber && gamma_cnt > wavenumber) break;//��alpha��gamma�Ĳ�����Ŀ�㹻ʱ����ֹѭ��
		t1->GetEntry(iwaveform);
		Float_t energy = Cal_Parameter[0]+Cal_Parameter[1]*TotalCharge+Cal_Parameter[2]*TotalCharge*TotalCharge;
		if (energy<700 && energy>620)//������620-700keV
		{
			gamma_flag = true;
		}
		else if (energy<2500 && energy>1800)//������1800-2500keV
		{
			gamma_flag = false;
		}
		else
		{
			continue;
		}
		if (gamma_flag)
		{
			if(gamma_cnt> = wavenumber) continue;
			gamma_cnt++;
			for (Int_t i=0;i<400;i++)
			{
				NormalizedTime[i]=i*0.4;
				NormalizedWave[i] = AlinedWave[i]/TotalCharge;
				TotalWave_gamma[i]+= NormalizedWave[i];
			}		
		}
		else
		{	
			if(alpha_cnt>= wavenumber) continue;
			alpha_cnt++;
			for (Int_t i=0;i<400;i++)
			{
				NormalizedTime[i]=i*0.4;
				NormalizedWave[i] = AlinedWave[i]/TotalCharge;
				TotalWave_alpha[i]+= NormalizedWave[i];
			}	
		}
		if(iwaveform%1000 == 0) cout<<"alpha_cnt="<<alpha_cnt<<" gamma_cnt="<<gamma_cnt<<endl;
	}
	cout<<"alpha_cnt="<<alpha_cnt<<" TotalWave_alpha[100]="<< TotalWave_alpha[100]<<" gamma_cnt="<<gamma_cnt<<" TotalWave_gamma[100]="<< TotalWave_gamma[100]<<endl;
	Float_t Average_alpha_temp[400];//ƽ������
	Float_t Average_gamma_temp[400];
	//���ı��ļ���ʽ��ƽ���������
	for (Int_t i=0;i<400;i++)
	{
		Average_alpha_temp[i] = TotalWave_alpha[i]/(Float_t)wavenumber;
		Average_gamma_temp[i] = TotalWave_gamma[i]/(Float_t)wavenumber;
		//cout << Average_alpha[i]<<"  "<< Average_gamma[i]<<endl;
	}
	//ƽ���˲���100��5��S-Kƽ��
	for (s=0;s<100;s++)
	{
		for (Int_t i=2;i<400-2;i++)
		{
			Float_t a[5];
			Float_t b[5];
			for (Int_t j=0;j<5;j++)
			{
				a[j]=Average_alpha_temp[i+j-2];
				b[j]=Average_gamma_temp[i+j-2];
			}
			Average_alpha[i] = (-3*a[0]+12*a[1]+17*a[2]+12*a[3]-3*a[4])/35;
			Average_gamma[i] = (-3*b[0]+12*b[1]+17*b[2]+12*b[3]-3*b[4])/35;
		}
		for (Int_t i=0;i<400;i++)
		{
			Average_alpha_temp[i] = Average_alpha[i];
			Average_gamma_temp[i] = Average_gamma[i];
		}
	}
	//���ı��ļ���ʽ��ƽ���������
	for (Int_t i=0;i<400;i++)
	{
		f_out_alpha <<setiosflags(ios_base::left)<<setw(10)<< i*0.4<<resetiosflags(ios_base::left) <<setw(10)<<Average_alpha[i]<<endl;
		f_out_gamma <<setiosflags(ios_base::left)<<setw(10)<< i*0.4<<resetiosflags(ios_base::left) <<setw(10)<<Average_gamma[i]<<endl;
		//cout << Average_alpha[i]<<"  "<< Average_gamma[i]<<endl;
	}
	//����ƽ��֮��ĵ��ͣ���һ������
	t2->Fill();
	t2->Write();

	gROOT->cd();//��������Ҫ������
	TCanvas *c1 = new TCanvas("c4","BC canvas",960,0,1000,600);
	TGraph *g_gamma = new TGraph(400,Time,Average_gamma);
	g_gamma->SetTitle("averaged waveforms;Time(ns);Amplitude(A.U.)");
	g_gamma->GetXaxis()->CenterTitle(kTRUE); 
	g_gamma->GetYaxis()->CenterTitle(kTRUE); 
	g_gamma->GetYaxis()->SetTitleOffset(1.2); 
	g_gamma->SetLineColor(kRed);
	g_gamma->SetLineWidth(1);
	g_gamma->Draw("AL");
	TGraph *g_alpha = new TGraph(400,Time,Average_alpha);
	g_alpha->SetLineColor(1);
	g_alpha->SetLineWidth(2);
	g_alpha->SetLineStyle(7);
	g_alpha->Draw("same");

	//��ͼ�л��Ʋο���ǣ��ҳ����߽��㣬�������ͼ��
	for (Int_t i=50;i<250;i++)
	{
		if((Average_alpha[i-1] > Average_gamma[i-1] && Average_alpha[i+1] < Average_gamma[i+1]) || (Average_alpha[i-1] < Average_gamma[i-1] && Average_alpha[i+1] > Average_gamma[i+1]))
		{
			char sf[100];
			sprintf(sf,"%.1f",Time[i]);  // float �� char
			TLatex text(Time[i]+1,Average_alpha[i],sf);
			text.SetTextSize(0.04);
			text.DrawClone();
			TLine* gLine = new TLine(Time[i],0,Time[i],0.02);
			gLine->SetLineStyle(7);
			gLine->Draw();
			i+=12;//�ҵ���������һЩ�㣬Ѱ����һ�����㣬��ֹ���ο���̫��
		}
	}
	//����legendͼ��
	TLegend *legend=new TLegend(0.75,0.78,0.90,0.88);//
	legend->SetTextFont(72);
	legend->SetTextSize(0.04);
	legend->AddEntry(g_alpha,"alpha","l");
	legend->AddEntry(g_gamma,"gamma","l");
	legend->Draw();
	c1->SaveAs(Form("%s/Picture/AveragedWaveforms.png",filepath));
	//���Ƶ��Ͳ��εĲ���ͼ
	Float_t alpha_max = TMath::MaxElement(400,Average_alpha);
	for (Int_t i=0;i<400;i++)
	{
		diff[i] =  100*(Average_alpha[i]-Average_gamma[i])/alpha_max;
	}
	TCanvas *c2 = new TCanvas("c2","Diff canvas",960,400,1000,600);
	TGraph *g_diff = new TGraph(400,Time,diff);
	g_diff->SetTitle(";Time(ns);Percent Difference(%)");
	g_diff->SetLineColor(1);
	g_diff->SetLineWidth(3);
	g_diff->GetXaxis()->CenterTitle(kTRUE); 
	g_diff->GetYaxis()->CenterTitle(kTRUE); 
	g_diff->GetYaxis()->SetTitleSize(0.05); 
	g_diff->GetYaxis()->SetTitleOffset(0.8);
	g_diff->GetXaxis()->SetTitleSize(0.05); 
	g_diff->GetXaxis()->SetTitleOffset(0.8);
	g_diff->Draw("AL");
	for (Int_t i=50;i<250;i++)
	{
		if(diff[i]*diff[i+1]<0 )
		{
			char sf[100];
			sprintf(sf,"%.1f",Time[i]);  // float �� char
			TLatex text1(Time[i]+1,Average_alpha[i]+0.05,sf);
			text1.SetTextSize(0.04);
			text1.DrawClone();
			TLine* gLine1 = new TLine(Time[i],-1.5,Time[i],2);
			gLine1->SetLineStyle(7);
			gLine1->Draw();
			i+=12;
		}
	}
	TF1 *f1 = new TF1("f1","0",0,160);	
	f1->SetLineColor(1);
	f1->SetLineStyle(9);//����
	f1->Draw("same");
	c2->SaveAs(Form("%s/Picture/Diff_average.png",filepath));
	f_out.Close();
	f_in.Close();
}
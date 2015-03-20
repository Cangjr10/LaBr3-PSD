//Updata time:2014/11/28
//Author�� cjr
//description: read out data from .txt files into .root files.
void ToTree()
{
    gROOT->Reset();    
/****����ṹ�����ڴ洢���ε�ԭʼ����************************************/
/****�ı��ļ���ʱ�䡢�������ݹ�1252��************************************/
    struct data_s
	{
		Float_t time[1252];
		Float_t ampl[1252];
	};
	data_s data_in;

    Int_t  iwaveform = 0;
    ifstream f_in;
	//���tree�ļ���Ŀ¼
    TFile *f_out = new TFile("C:/root/Projects/LaBr3/SGB_2inch_background_WithShield/Tree/Data0.root","recreate");
    TTree *t0 = new TTree("t0","data");
    t0->Branch("data",&data_in.time,"time[1252]/F:ampl[1252]");

    Char_t  title[50];
    Float_t time_t,ampl_t;
    Int_t   iline;
	//���ε�����Ŀ���ļ��������������
	const int txt_number = 99999;
	for(iwaveform=0;iwaveform<txt_number+1;iwaveform++)
	{
		//text�����ļ���Ŀ¼
		f_in.open(Form("F:/Nuclear Electronics/LaBr3 PSD/LaBr3WaveData/SGB_2inch_shield_background/LaBr3_20141113_HighEnergy1/C12inch_LaBr_20141113_%05d.txt",iwaveform),ios_base::in); 
		//����ǰ5�е��ļ�˵��
		f_in.getline(title,50);
		f_in.getline(title,50);
		f_in.getline(title,50);
		f_in.getline(title,50);
		f_in.getline(title,50);
        iline = 0;
		while(1)
		{
			f_in >> time_t >> ampl_t;
			if (f_in.get()==EOF) break;//�ļ���β
			data_in.time[iline] = time_t*1000000000; //ns
			data_in.ampl[iline] = ampl_t*1000;       //mV
		    iline++;
		}
        t0->Fill();
		f_in.close();
		f_in.clear();
		if(iwaveform%100==0) cout <<  iwaveform << " files" << " processed." << endl;
	}
	t0->Write();
	f_out->Close();
	cout << "Tree file generated." <<endl;
}
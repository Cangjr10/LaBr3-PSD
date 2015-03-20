#include <cmath>
void Correlation_CCM(char* filepath = "C:/ROOT/Projects/LaBr3/CFD/SGB_2inch_background_WithShield")
{
	TH2D *h2D_alpha= new TH2D("h2D_alpha","TotalCharge VS PartialCharge;TotalCharge;PartialCharge",100,0,9000,1000,0,6000);
	TH2D *h2D_gamma= new TH2D("h2D_gamma","TotalCharge VS PartialCharge;TotalCharge;PartialCharge",100,800,500,100,400,3000);
	TH1D *h1D_CCM_gamma = new TH1D("h1D_gamma","CCM Distribution",50,0.635,0.665);
	//Energy Calibration
	Double_t Cal_Parameter[3];
	TFile f_in1(Form("%s/Tree/EnergyCalibration.root",filepath));
	TTree *t_energy_cali = (TTree*) f_in1.GetObjectChecked("t_energy_cali","TTree");
	TBranch *b_cali = t_energy_cali->GetBranch("Cal_Parameter");
	b_cali->SetAddress(Cal_Parameter);
	t_energy_cali->GetEntry(0);
	//CCM info ,Partial and Total charge.
	TFile f_in(Form("%s/Tree/CCM.root",filepath));
	TTree *t_CCM = (TTree*) f_in.GetObjectChecked("t_CCM","TTree");
	Long64_t nentries = t_CCM->GetEntries();
	Float_t Q_total,Q_part;
	t_CCM->SetBranchAddress("Q_part",&Q_part);
	t_CCM->SetBranchAddress("Q_total",&Q_total);
	Long64_t nEntries = t_CCM->GetEntries();
	//h2D_gamma->Sumw2(0);
	Int_t Count = 0;
	for(Int_t iEntry=0;iEntry < nEntries;iEntry++)
	{
		t_CCM->GetEntry(iEntry);
		Float_t Energy0 = Cal_Parameter[0]+Cal_Parameter[1]*Q_total+Cal_Parameter[2]*Q_total*Q_total;
		if (Energy0 < 500 && Energy0 >200)
		{
			Count++;
			if(Count>1000) break;
			h2D_gamma->Fill(Q_total,Q_part);
			h1D_CCM_gamma->Fill(Q_part/Q_total);
		}
		else if (Energy0>1600 && Energy0<2500)
		{
			h2D_alpha->Fill(Q_total,Q_part);
		}
	}
 	TH1D * h_Qp = h2D_gamma->ProjectionY();
 	TH1D * h_Qt = h2D_gamma->ProjectionX();
	Double_t sigma_Qp = h_Qp->GetStdDev();
	Double_t sigma_Qt = h_Qt->GetStdDev();
	Double_t Mean_Qp = h_Qp->GetMean();
	Double_t Mean_Qt = h_Qt->GetMean();
	cout<<"Mean of Qp is: "<<Mean_Qp<<endl;
	cout<<"Mean of Qt is: "<<Mean_Qt<<endl;
	cout<<"standard deviation of Qp is: "<<sigma_Qp<<endl;
	cout<<"standard deviation of Qt is: "<<sigma_Qt<<endl;
	Double_t Covariance = h2D_gamma->GetCovariance();
	cout<<"Covariance of Gamma Event is : "<<Covariance<<endl;
	Double_t Sigma_CCM = TMath::Sqrt(pow(sigma_Qp,2)/pow(Mean_Qt,2)+pow(Mean_Qp,2)*pow(sigma_Qt,2)/pow(Mean_Qt,4)-2*Mean_Qp*Covariance/pow(Mean_Qt,3));
	cout<<"uncertainty of CCM is "<<Sigma_CCM<<endl;
	Double_t CorrelationFactor = h2D_gamma->GetCorrelationFactor();
	cout<<"Correlation Factor of Gamma Event is : "<<CorrelationFactor<<endl;
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1","A canvas",0,0,700,400);
	h2D_gamma->Draw("");
	TCanvas *c2 = new TCanvas("c2","B canvas",0,500,700,400);
	h1D_CCM_gamma->Draw();
	h1D_CCM_gamma->Fit("gaus");
	gROOT->cd();
	TCanvas *c3 = new TCanvas("c3","C canvas",750,0,700,400);
// 	TH1D * projh2X = h2D_gamma->ProjectionY();
// 	projh2X->Draw();
	h2D_gamma->ProjectionX()->Draw();
	TCanvas *c4 = new TCanvas("c4","D canvas",750,500,700,400);
	// 	TH1D * projh2X = h2D_gamma->ProjectionY();
	// 	projh2X->Draw();
	h2D_gamma->ProjectionY()->Draw();
// 	//  		c1->Divide(1,2);
// 	//  		c1->cd(1);
// 	// c1->SetLogy();
// 	// c1->SetLogx();
// 	TF1 f("Linear law","x*[1]",0,4000);//[0]+
// 	// Let's make the function line nicer
// 	f.SetLineColor(kBlue); f.SetLineStyle(2);
// 	// Fit it to the graph and draw it
// 	h2D_gamma->Fit(&f);
// 	TF1 f1("Linear law","x*[1]",4000,9000);//[0]+
// 	// Let's make the function line nicer
// 	f1.SetLineColor(kRed); f1.SetLineStyle(2);
// 	// Fit it to the graph and draw it
// 	h2D_alpha->Fit(&f1);
// 	h2D_gamma->SetMarkerColor(kBlue);
// 	h2D_gamma->SetMarkerSize(0.3);
// 	h2D_gamma->SetMarkerStyle(21);
// 	h2D_gamma->Draw();
// 	f.DrawClone("ACsame");
// 	f1.SetLineWidth(3);
// 	f1.Draw("same");
// 	h2D_alpha->SetMarkerColor(kRed);
// 	h2D_alpha->SetMarkerSize(0.3);
// 	h2D_alpha->SetMarkerStyle(21);
// 	h2D_alpha->Draw("same");
// 
// 	TLegend *legend1=new TLegend(0.75,0.78,0.90,0.88);//
// 	legend1->SetTextFont(72);
// 	legend1->SetTextSize(0.04);
// 	legend1->AddEntry(h2D_alpha,"alpha_red","P");
// 	legend1->AddEntry(h2D_gamma,"gamma_blue","P");
// 	legend1->Draw();
}
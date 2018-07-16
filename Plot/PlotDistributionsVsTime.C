#include "PROSPECT_Style.cc"
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPave.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TColor.h"
#include "TExec.h"
#include <sstream>
#include "TLatex.h"
#include "TMath.h"
#include "TGraphErrors.h"

int PlotDistributionsVsTime(int timeBin){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	TFile *fg = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fg){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	
	TGraphErrors *grRate    = (TGraphErrors*)fg->Get("grRate");
	fg->Close();

	//-------------------------------------------------------------------------------------------------------
	TFile *f = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!f){
		printf("File not found. Exiting. \n");
		return -1;
	}
cout<<"Loading hists"<<endl;
	TH1F *hRnPoDt 	 	= (TH1F*)f->Get(Form("hRnPoDt_%i",timeBin));
	TF1  *fRnPoDtExp 	= hRnPoDt->GetFunction("fRnPoDtExp");
	fRnPoDtExp->SetRange(0,2.57*5.0);

	TH1F *hRnPSD 	 	= (TH1F*)f->Get(Form("hRnPSD_%i",timeBin));
	TF1  *fRnPSDGaus 	= hRnPSD->GetFunction("fRnPSDGaus");

	TH1F *hPoPSD 	 	= (TH1F*)f->Get(Form("hPoPSD_%i",timeBin));
	TF1  *fPoPSDGaus 	= hPoPSD->GetFunction("fPoPSDGaus");

	TH1F *hRnEn   		= (TH1F*)f->Get(Form("hRnEn_%i",timeBin));
	TF1  *fRnEnCB 		= hRnEn->GetFunction("fRnEnCB");

	TH1F *hPoEn 		= (TH1F*)f->Get(Form("hPoEn_%i",timeBin));
	TF1  *fPoEnGaus 	= hPoEn->GetFunction("fPoEnGaus");

	TH1F *hRnPos 		= (TH1F*)f->Get(Form("hRnPos_%i",timeBin));

	TH1F *hPoPos 		= (TH1F*)f->Get(Form("hPoPos_%i",timeBin));

	TH1F *hRnPoDz 	  	= (TH1F*)f->Get(Form("hRnPoDz_%i",timeBin));
	TF1  *fRnPoDzGaus 	= hRnPoDz->GetFunction("fRnPoDzGaus");

	TH2F *hRnPoPSDvsEn 	= (TH2F*)f->Get(Form("hRnPoPSDvsEn_%i",timeBin));

	TH2F *hPoEnvsRnEn 	= (TH2F*)f->Get(Form("hPoEnvsRnEn_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	int p_col = 8, d_col = 9;

cout<<"Plotting hists"<<endl;
	TCanvas *cRnPoDt = new TCanvas("cRnPoDt","RnPo Dt",1);
	gPad->SetGrid();
	hRnPoDt->SetMarkerSize(0.5);
	hRnPoDt->Draw();
	hRnPoDt->GetXaxis()->SetTitle("t_{#alphaPo} - t_{#alphaRn} [ms]");
	fRnPoDtExp->SetLineStyle(2);
	fRnPoDtExp->Draw("same");
	pt = new TPaveText(0.68,0.63,0.95,0.9,"NDCNB");
	pt->AddText(Form("Entries    %.0f",hRnPoDt->GetEntries()));
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDtExp->GetChisquare(),fRnPoDtExp->GetNDF()));
	pt->AddText(Form("Prob    %f",fRnPoDtExp->GetProb()));
	pt->AddText(Form("N    %.2f #pm %.2f",fRnPoDtExp->GetParameter(0),fRnPoDtExp->GetParError(0)));
	pt->AddText(Form("#tau_{Po^{215}}    %.2f #pm %.2f ms",fRnPoDtExp->GetParameter(1),fRnPoDtExp->GetParError(1)));
	pt->AddText(Form("R    %.2f #pm %.2f Hz",grRate->GetY()[timeBin],grRate->GetEY()[timeBin]));
	pt->Draw();
	cRnPoDt->SaveAs(Form("%s/RnPoDt_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoDt->SaveAs(Form("%s/RnPoDt_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));


	TCanvas *cRnPoPSD = new TCanvas("cRnPoPSD","RnPo PSD",1);
	gPad->SetGrid();
	hRnPSD->SetLineColor(p_col);
	hPoPSD->SetLineColor(d_col);
	hRnPSD->SetMarkerColor(p_col);
	hPoPSD->SetMarkerColor(d_col);
	hPoPSD->Draw();		
	hRnPSD->Draw("same");
	fRnPSDGaus->SetLineStyle(2);
	fPoPSDGaus->SetLineStyle(2);
	fRnPSDGaus->Draw("same");
	fPoPSDGaus->Draw("same");
	pt = new TPaveText(0.68,0.55,0.95,0.9,"NDCNB");
	TText *tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPSDGaus->GetChisquare(),fRnPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.3f",fRnPSDGaus->GetProb()));
	pt->AddText(Form("#mu    %.2f #pm %.2f",fRnPSDGaus->GetParameter(1),fRnPSDGaus->GetParError(1)));
	pt->AddText(Form("#sigma    %.2f #pm %.2f",fRnPSDGaus->GetParameter(2),fRnPSDGaus->GetParError(2)));
	pt->AddText("");
	TText *tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoPSDGaus->GetChisquare(),fPoPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.3f",fPoPSDGaus->GetProb()));
	pt->AddText(Form("#mu    %.2f #pm %.2f",fPoPSDGaus->GetParameter(1),fPoPSDGaus->GetParError(1)));
	pt->AddText(Form("#sigma    %.2f #pm %.2f",fPoPSDGaus->GetParameter(2),fPoPSDGaus->GetParError(2)));
	tRn->SetTextColor(p_col);
	tPo->SetTextColor(d_col);
	pt->Draw();
	cRnPoPSD->SaveAs(Form("%s/RnPoPSD_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoPSD->SaveAs(Form("%s/RnPoPSD_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	
	TCanvas *cRnPoEn = new TCanvas("cRnPoEn","RnPo En",1);
	gPad->SetGrid();
	hPoEn->GetYaxis()->SetTitleOffset(0.9);
	hRnEn->SetLineColor(p_col);
	hPoEn->SetLineColor(d_col);
	hRnEn->SetMarkerColor(p_col);
	hPoEn->SetMarkerColor(d_col);
	hPoEn->Draw();
	hRnEn->Draw("same");
	fRnEnCB->SetLineStyle(2);
	fPoEnGaus->SetLineStyle(2);
	fRnEnCB->Draw("same");
	fPoEnGaus->Draw("same");
	pt = new TPaveText(0.68,0.55,0.95,0.9,"NDCNB");
	tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnEnCB->GetChisquare(),fRnEnCB->GetNDF()));
	pt->AddText(Form("Prob    %.3f",fRnEnCB->GetProb()));
	pt->AddText(Form("#mu    %.3f #pm %.3f",fRnEnCB->GetParameter(1),fRnEnCB->GetParError(1)));
	pt->AddText(Form("#sigma    %.3f #pm %.3f",fRnEnCB->GetParameter(2),fRnEnCB->GetParError(2)));
	pt->AddText("");
	tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoEnGaus->GetChisquare(),fPoEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.3f",fPoEnGaus->GetProb()));
	pt->AddText(Form("#mu    %.3f #pm %.3f",fPoEnGaus->GetParameter(1),fPoEnGaus->GetParError(1)));
	pt->AddText(Form("#sigma    %.3f #pm %.3f",fPoEnGaus->GetParameter(2),fPoEnGaus->GetParError(2)));
	tRn->SetTextColor(p_col);
	tPo->SetTextColor(d_col);
	pt->Draw();
	cRnPoEn->SaveAs(Form("%s/RnPoEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoEn->SaveAs(Form("%s/RnPoEn_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
		

	TCanvas *cRnPoDz = new TCanvas("cRnPoDz","RnPo Dz",1);
	gPad->SetGrid();
	hRnPoDz->GetYaxis()->SetTitleOffset(0.9);
	hRnPoDz->Draw();
	hRnPoDz->GetXaxis()->SetTitle("z_{#alphaPo} - z_{#alphaRn} [ms]");
	fRnPoDzGaus->SetLineStyle(2);
	fRnPoDzGaus->Draw("same");	
	pt = new TPaveText(0.68,0.63,0.95,0.9,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDzGaus->GetChisquare(),fRnPoDzGaus->GetNDF()));
	pt->AddText(Form("Prob    %.3f",fRnPoDzGaus->GetProb()));
	pt->AddText(Form("#mu    %.2f #pm %.2f",fRnPoDzGaus->GetParameter(1),fRnPoDzGaus->GetParError(1)));
	pt->AddText(Form("#sigma    %.2f #pm %.2f",fRnPoDzGaus->GetParameter(2),fRnPoDzGaus->GetParError(2)));
	pt->Draw();
	cRnPoDz->SaveAs(Form("%s/RnPoDz_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoDz->SaveAs(Form("%s/RnPoDz_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));

	
	TCanvas *cRnPoPos = new TCanvas("cRnPoPos","RnPo Position",1);
	gPad->SetGrid();
	hRnPos->SetLineColor(p_col);
	hPoPos->SetLineColor(d_col);
	hRnPos->GetYaxis()->SetTitleOffset(1.0);
	hRnPos->Draw();
	hPoPos->Draw("same");
	pt = new TPaveText(0.39,0.2,0.66,0.45,"NDCNB");
	tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#mu    %.2f #pm %.2f",hRnPos->GetMean(),hRnPos->GetMeanError()));
	pt->AddText(Form("RMS    %.2f #pm %.2f",hRnPos->GetRMS(),hRnPos->GetRMSError()));
	pt->AddText("");
	tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#mu    %.2f #pm %.2f",hPoPos->GetMean(),hPoPos->GetMeanError()));
	pt->AddText(Form("RMS    %.2f #pm %.2f",hPoPos->GetRMS(),hPoPos->GetRMSError()));
	tRn->SetTextColor(p_col);
	tPo->SetTextColor(d_col);
	pt->Draw();
	cRnPoPos->SaveAs(Form("%s/RnPoPos_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoPos->SaveAs(Form("%s/RnPoPos_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));


	TCanvas *cRnPoPSDvsEn = new TCanvas("cRnPoPSDvsEn","RnPo PSD vs En",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsEn->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsEn->Draw("colz");
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	
	TCanvas *cPoEnvsRnEn = new TCanvas("cPoEnvsRnEn","Po En vs Rn En",650,600);
	gPad->SetRightMargin(0.12);
	gPad->SetLeftMargin(0.12);
	hPoEnvsRnEn->GetYaxis()->SetTitleOffset(1.1);
	hPoEnvsRnEn->Draw("colz");
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_TimeBin%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_TimeBin%i.png",gSystem->Getenv("AD_AC227_PLOTS"),timeBin));


	return 0;
} 	//end PlotDistributionsVsTime

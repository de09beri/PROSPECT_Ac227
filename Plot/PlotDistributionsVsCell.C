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

int PlotDistributionsVsCell(int cellNum){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	const int NUMEXCLUDECELLS = 31;
	int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};	

	if(find(begin(ExcludeCellArr), end(ExcludeCellArr), cellNum) != end(ExcludeCellArr)){
		printf("No information for cell %i. Exiting. \n",cellNum);
		return -1;
	}

	int numLower = 0.0;
	for(int i=0;i<NUMEXCLUDECELLS;i++){
		int v = ExcludeCellArr[i];
		if(v<cellNum) numLower++;
		if(v>cellNum) break;
	}	

	int grIDX = cellNum - numLower;

	TFile *fg = new TFile(Form("%s/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fg){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	TGraphErrors *grRnPSDEff    = (TGraphErrors*)fg->Get("grRnPSDEff");
    TGraphErrors *grPoPSDEff    = (TGraphErrors*)fg->Get("grPoPSDEff");
    TGraphErrors *grRnEnEff     = (TGraphErrors*)fg->Get("grRnEnEff");
    TGraphErrors *grPoEnEff     = (TGraphErrors*)fg->Get("grPoEnEff");
    TGraphErrors *grRnPoDzEff   = (TGraphErrors*)fg->Get("grRnPoDzEff");
	fg->Close();

	cout<<grIDX<<endl;
	cout<<grRnPSDEff->GetX()[grIDX]<<endl;
	//-------------------------------------------------------------------------------------------------------
	TFile *f = new TFile(Form("%s/Ac227_HistsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!f){
		printf("File not found. Exiting. \n");
		return -1;
	}


	TH1F *hSelectDt = (TH1F*)f->Get(Form("hSelectDt_%i",cellNum));
	TH1F *hBGDt 	= (TH1F*)f->Get(Form("hBGDt_%i",cellNum));
	TH1F *hRnPoDt 	= (TH1F*)f->Get(Form("hRnPoDt_%i",cellNum));

	TF1  *fRnPoDtExp = hRnPoDt->GetFunction("fRnPoDtExp");
	fRnPoDtExp->SetRange(0,2.57*5.0);

	TH1F *hSelectPromptPSD = (TH1F*)f->Get(Form("hSelectPromptPSD_%i",cellNum));
	TH1F *hBGPromptPSD 	   = (TH1F*)f->Get(Form("hBGPromptPSD_%i",cellNum));
	TH1F *hRnPSD 	 	   = (TH1F*)f->Get(Form("hRnPSD_%i",cellNum));

	TF1  *fRnPSDGaus = hRnPSD->GetFunction("fRnPSDGaus");
	fRnPSDGaus->SetRange(0.15,0.37);

	TH1F *hSelectDelayPSD = (TH1F*)f->Get(Form("hSelectDelayPSD_%i",cellNum));
	TH1F *hBGDelayPSD 	  = (TH1F*)f->Get(Form("hBGDelayPSD_%i",cellNum));
	TH1F *hPoPSD 		  = (TH1F*)f->Get(Form("hPoPSD_%i",cellNum));

	TF1  *fPoPSDGaus = hPoPSD->GetFunction("fPoPSDGaus");
	fPoPSDGaus->SetRange(0.15,0.37);

	TH1F *hSelectPromptEn = (TH1F*)f->Get(Form("hSelectPromptEn_%i",cellNum));
	TH1F *hBGPromptEn 	  = (TH1F*)f->Get(Form("hBGPromptEn_%i",cellNum));
	TH1F *hRnEn 		  = (TH1F*)f->Get(Form("hRnEn_%i",cellNum));

	TF1  *fRnEnCB = hRnEn->GetFunction("fRnEnCB");
	fRnEnCB->SetRange(0.49,1.16);

	TH1F *hSelectDelayEn = (TH1F*)f->Get(Form("hSelectDelayEn_%i",cellNum));
	TH1F *hBGDelayEn 	 = (TH1F*)f->Get(Form("hBGDelayEn_%i",cellNum));
	TH1F *hPoEn 		 = (TH1F*)f->Get(Form("hPoEn_%i",cellNum));

	TF1  *fPoEnGaus = hPoEn->GetFunction("fPoEnGaus");
	fPoEnGaus->SetRange(0.49,1.16);

	TH1F *hSelectPromptPos = (TH1F*)f->Get(Form("hSelectPromptPos_%i",cellNum));
	TH1F *hBGPromptPos 	   = (TH1F*)f->Get(Form("hBGPromptPos_%i",cellNum));
	TH1F *hRnPos 		   = (TH1F*)f->Get(Form("hRnPos_%i",cellNum));

	TH1F *hSelectDelayPos = (TH1F*)f->Get(Form("hSelectDelayPos_%i",cellNum));
	TH1F *hBGDelayPos 	  = (TH1F*)f->Get(Form("hBGDelayPos_%i",cellNum));
	TH1F *hPoPos 		  = (TH1F*)f->Get(Form("hPoPos_%i",cellNum));

	TH1F *hSelectDz = (TH1F*)f->Get(Form("hSelectDz_%i",cellNum));
	TH1F *hBGDz 	= (TH1F*)f->Get(Form("hBGDz_%i",cellNum));
	TH1F *hRnPoDz 	= (TH1F*)f->Get(Form("hRnPoDz_%i",cellNum));

	TF1  *fRnPoDzGaus = hRnPoDz->GetFunction("fRnPoDzGaus");

	TH2F *hSelectPSDvsEn = (TH2F*)f->Get(Form("hSelectPSDvsEn_%i",cellNum));
	TH2F *hBGPSDvsEn 	 = (TH2F*)f->Get(Form("hBGPSDvsEn_%i",cellNum));
	TH2F *hRnPoPSDvsEn 	 = (TH2F*)f->Get(Form("hRnPoPSDvsEn_%i",cellNum));

	TH2F *hSelectPSDvsPos = (TH2F*)f->Get(Form("hSelectPSDvsPos_%i",cellNum));
	TH2F *hBGPSDvsPos 	  = (TH2F*)f->Get(Form("hBGPSDvsPos_%i",cellNum));
	TH2F *hRnPoPSDvsPos   = (TH2F*)f->Get(Form("hRnPoPSDvsPos_%i",cellNum));

	TH2F *hSelectEnvsPos = (TH2F*)f->Get(Form("hSelectEnvsPos_%i",cellNum));
	TH2F *hBGEnvsPos 	 = (TH2F*)f->Get(Form("hBGEnvsPos_%i",cellNum));
	TH2F *hRnPoEnvsPos   = (TH2F*)f->Get(Form("hRnPoEnvsPos_%i",cellNum));

	TH2F *hSelectDelayEnvsPromptEn = (TH2F*)f->Get(Form("hSelectDelayEnvsPromptEn_%i",cellNum));
	TH2F *hBGDelayEnvsPromptEn 	   = (TH2F*)f->Get(Form("hBGDelayEnvsPromptEn_%i",cellNum));
	TH2F *hPoEnvsRnEn 			   = (TH2F*)f->Get(Form("hPoEnvsRnEn_%i",cellNum));

	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	//int p_col = 8, d_col = 9;
	int p_col = 1, d_col = 4;

	TCanvas *cRnPoDt = new TCanvas("cRnPoDt","RnPo Dt",1);
	gPad->SetGrid();
	gPad->SetLogy();
	hSelectDt->SetLineColor(8);
	hSelectDt->SetLineWidth(1);
	hSelectDt->SetMinimum(0.1);
	hSelectDt->Draw("HIST");
	hBGDt->SetLineColor(kBlack);
	hBGDt->SetLineWidth(1);
	hBGDt->Draw("same&HIST");
	hRnPoDt->SetLineWidth(1);
	hRnPoDt->Draw("same");
	hRnPoDt->GetXaxis()->SetTitle("t_{#alphaPo} - t_{#alphaRn} [ms]");
	fRnPoDtExp->SetLineStyle(2);
	fRnPoDtExp->SetRange(0,12.85);
	fRnPoDtExp->Draw("same");
	pt = new TPaveText(0.68,0.63,0.95,0.9,"NDCNB");
	pt->AddText(Form("Entries    %.0f",hRnPoDt->GetEntries()));
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDtExp->GetChisquare(),fRnPoDtExp->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDtExp->GetProb()));
	pt->AddText(Form("N    %.2f #pm %.2f",fRnPoDtExp->GetParameter(0),fRnPoDtExp->GetParError(0)));
	pt->AddText(Form("#tau_{Po^{215}}    %.2f #pm %.2f ms",fRnPoDtExp->GetParameter(1),fRnPoDtExp->GetParError(1)));
	pt->Draw();
	cRnPoDt->SaveAs(Form("%s/RnPoDt_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoDt->SaveAs(Form("%s/RnPoDt_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));


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
	pt = new TPaveText(0.68,0.59,0.95,0.9,"NDCNB");
	TText *tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPSDGaus->GetChisquare(),fRnPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grRnPSDEff->GetY()[grIDX],grRnPSDEff->GetEY()[grIDX]));
	pt->AddText("");
	TText *tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoPSDGaus->GetChisquare(),fPoPSDGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoPSDGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grPoPSDEff->GetY()[grIDX],grPoPSDEff->GetEY()[grIDX]));
	tRn->SetTextColor(p_col);
	tPo->SetTextColor(d_col);
	pt->Draw();
	cRnPoPSD->SaveAs(Form("%s/RnPoPSD_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPSD->SaveAs(Form("%s/RnPoPSD_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	
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
	pt = new TPaveText(0.68,0.59,0.95,0.9,"NDCNB");
	tRn = pt->AddText("Rn^{219}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnEnCB->GetChisquare(),fRnEnCB->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnEnCB->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grRnEnEff->GetY()[grIDX],grRnEnEff->GetEY()[grIDX]));
	pt->AddText("");
	tPo = pt->AddText("Po^{215}");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fPoEnGaus->GetChisquare(),fPoEnGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fPoEnGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grPoEnEff->GetY()[grIDX],grPoEnEff->GetEY()[grIDX]));
	tRn->SetTextColor(p_col);
	tPo->SetTextColor(d_col);
	pt->Draw();
	cRnPoEn->SaveAs(Form("%s/RnPoEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoEn->SaveAs(Form("%s/RnPoEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
		

	TCanvas *cRnPoDz = new TCanvas("cRnPoDz","RnPo Dz",1);
	gPad->SetGrid();
	hRnPoDz->GetYaxis()->SetTitleOffset(0.9);
	hRnPoDz->Draw();
	hRnPoDz->GetXaxis()->SetTitle("z_{#alphaPo} - z_{#alphaRn} [ms]");
	fRnPoDzGaus->SetLineStyle(2);
	fRnPoDzGaus->Draw("same");	
	pt = new TPaveText(0.68,0.63,0.95,0.9,"NDCNB");
	pt->AddText(Form("#Chi^{2}/NDF    %.2f/%d",fRnPoDzGaus->GetChisquare(),fRnPoDzGaus->GetNDF()));
	pt->AddText(Form("Prob    %.2f",fRnPoDzGaus->GetProb()));
	pt->AddText(Form("Eff    %.2f #pm %.2f",grRnPoDzEff->GetY()[grIDX],grRnPoDzEff->GetEY()[grIDX]));
	pt->Draw();
	cRnPoDz->SaveAs(Form("%s/RnPoDz_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoDz->SaveAs(Form("%s/RnPoDz_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));


	TCanvas *cRnPoPos = new TCanvas("cRnPoPos","RnPo Position",1);
    gPad->SetGrid();
    hRnPos->SetLineColor(p_col);
    hRnPos->SetMarkerColor(p_col);
    hPoPos->SetLineColor(d_col);
    hPoPos->SetMarkerColor(d_col);
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
	cRnPoPos->SaveAs(Form("%s/RnPoPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPos->SaveAs(Form("%s/RnPoPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));


	TCanvas *cRnPoPSDvsEn = new TCanvas("cRnPoPSDvsEn","RnPo PSD vs En",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsEn->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsEn->Draw("colz");
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	
	TCanvas *cRnPoPSDvsPos = new TCanvas("cRnPoPSDvsPos","RnPo PSD vs Pos",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoPSDvsPos->GetYaxis()->SetTitleOffset(1.1);
	hRnPoPSDvsPos->Draw("colz");
	cRnPoPSDvsPos->SaveAs(Form("%s/RnPoPSDvsPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoPSDvsPos->SaveAs(Form("%s/RnPoPSDvsPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	TCanvas *cRnPoEnvsPos = new TCanvas("cRnPoEnvsPos","RnPo En vs Pos",1);
	gPad->SetRightMargin(0.1);
	gPad->SetLeftMargin(0.12);
	hRnPoEnvsPos->GetYaxis()->SetTitleOffset(1.1);
	hRnPoEnvsPos->GetYaxis()->SetTitle("Energy [MeVee]");
	hRnPoEnvsPos->Draw("colz");
	cRnPoEnvsPos->SaveAs(Form("%s/RnPoEnvsPos_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cRnPoEnvsPos->SaveAs(Form("%s/RnPoEnvsPos_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));

	TCanvas *cPoEnvsRnEn = new TCanvas("cPoEnvsRnEn","Po En vs Rn En",650,600);
	gPad->SetRightMargin(0.12);
	gPad->SetLeftMargin(0.12);
	hPoEnvsRnEn->GetYaxis()->SetTitleOffset(1.1);
	hPoEnvsRnEn->Draw("colz");
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_Cell%i.pdf",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));
	cPoEnvsRnEn->SaveAs(Form("%s/PoEnvsRnEn_Cell%i.png",gSystem->Getenv("AD_AC227_PLOTS"),cellNum));


	return 0;
} 	//end PlotDistributionsVsCell

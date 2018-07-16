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

TH1F *makeHist(TGraphErrors *gr){
	int numPt = gr->GetN();
	double grx, gry, grxErr, gryErr;

	TH1F *h = new TH1F("h","h",154,0,154);
	
	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		h->SetBinContent((int)grx,gry);
		h->SetBinError((int)grx,gryErr);
	}

	return h;
}

TH1F *makeEnHist(TGraphErrors *gr){
	int numPt = gr->GetN();
	double grx, gry, grxErr, gryErr;

	TH1F *h = new TH1F("h","h",154,0,154);
	
	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		gry = gry*1000.0;		//convert MeV to keV
		gryErr = gryErr*1000.0;

		h->SetBinContent((int)grx,gry);
		h->SetBinError((int)grx,gryErr);
	}

	return h;
}

int PlotRxOnVsRxOff(){
	int timeBin = 0;

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	//-------------------------------------------------------------------------------------------------------
	TFile *fOn = new TFile(Form("%s/RxOn/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!fOn){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_On 	= (TH1F*)fOn->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_On 	= (TH1F*)fOn->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_On   	= (TH1F*)fOn->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_On 		= (TH1F*)fOn->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_On 	= (TH1F*)fOn->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_On 	= (TH1F*)fOn->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_On 	= (TH1F*)fOn->Get(Form("hRnPoDz_%i",timeBin));


	TFile *fOff = new TFile(Form("%s/RxOff/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!fOff){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TH1F *hRnPSD_Off 	= (TH1F*)fOff->Get(Form("hRnPSD_%i",timeBin));
	TH1F *hPoPSD_Off 	= (TH1F*)fOff->Get(Form("hPoPSD_%i",timeBin));
	TH1F *hRnEn_Off   	= (TH1F*)fOff->Get(Form("hRnEn_%i",timeBin));
	TH1F *hPoEn_Off 	= (TH1F*)fOff->Get(Form("hPoEn_%i",timeBin));
	TH1F *hRnPos_Off 	= (TH1F*)fOff->Get(Form("hRnPos_%i",timeBin));
	TH1F *hPoPos_Off 	= (TH1F*)fOff->Get(Form("hPoPos_%i",timeBin));
	TH1F *hRnPoDz_Off 	= (TH1F*)fOff->Get(Form("hRnPoDz_%i",timeBin));

	//-------------------------------------------------------------------------------------------------------
	TFile *fgrOn = new TFile(Form("%s/RxOn/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!fgrOn){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grRate_On 		= (TGraphErrors*)fgrOn->Get("grRate");
	TGraphErrors *grPoEnMean_On 	= (TGraphErrors*)fgrOn->Get("grPoEnMean");
	TGraphErrors *grPoEnSigma_On 	= (TGraphErrors*)fgrOn->Get("grPoEnSigma"); 
	TGraphErrors *grRnPoDzMean_On 	= (TGraphErrors*)fgrOn->Get("grRnPoDzMean");
	TGraphErrors *grRnPoDzSigma_On 	= (TGraphErrors*)fgrOn->Get("grRnPoDzSigma");

	TH1F *hRate_On 		  = makeHist(grRate_On); 
	TH1F *hPoEnMean_On    = makeEnHist(grPoEnMean_On);
	TH1F *hPoEnSigma_On   = makeEnHist(grPoEnSigma_On);
	TH1F *hRnPoDzMean_On  = makeHist(grRnPoDzMean_On);
	TH1F *hRnPoDzSigma_On = makeHist(grRnPoDzSigma_On);

	TFile *fgrOff = new TFile(Form("%s/RxOff/Ac227_GraphsPerCell.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));
	if(!fgrOff){
		printf("File not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grRate_Off 		= (TGraphErrors*)fgrOff->Get("grRate");
	TGraphErrors *grPoEnMean_Off 	= (TGraphErrors*)fgrOff->Get("grPoEnMean");
	TGraphErrors *grPoEnSigma_Off 	= (TGraphErrors*)fgrOff->Get("grPoEnSigma"); 
	TGraphErrors *grRnPoDzMean_Off 	= (TGraphErrors*)fgrOff->Get("grRnPoDzMean");
	TGraphErrors *grRnPoDzSigma_Off = (TGraphErrors*)fgrOff->Get("grRnPoDzSigma");

	TH1F *hRate_Off 	   = makeHist(grRate_Off);
	TH1F *hPoEnMean_Off    = makeEnHist(grPoEnMean_Off);
	TH1F *hPoEnSigma_Off   = makeEnHist(grPoEnSigma_Off);
	TH1F *hRnPoDzMean_Off  = makeHist(grRnPoDzMean_Off);
	TH1F *hRnPoDzSigma_Off = makeHist(grRnPoDzSigma_Off);


	double RxOnTime = 792.784039;	//hrs
	double RxOffTime = 655.099961;	//hrs

	//-------------------------------------------------------------------------------------------------------
	hRnEn_On->Scale(1/RxOnTime);
	hPoEn_On->Scale(1/RxOnTime);
	hRnPoDz_On->Scale(1/RxOnTime);
	
	hRnEn_Off->Scale(1/RxOffTime);
	hPoEn_Off->Scale(1/RxOffTime);
	hRnPoDz_Off->Scale(1/RxOffTime);

	//-------------------------------------------------------------------------------------------------------
	TH1F *hRnEn_Sub = (TH1F*)hRnEn_On->Clone();
	hRnEn_Sub->Sumw2();
	hRnEn_Sub->Add(hRnEn_Off,-1);
	hRnEn_Sub->Divide(hRnEn_Off);

	TH1F *hPoEn_Sub = (TH1F*)hPoEn_On->Clone();
	hPoEn_Sub->Sumw2();
	hPoEn_Sub->Add(hPoEn_Off,-1);
	hPoEn_Sub->Divide(hPoEn_Off);

	TH1F *hRnPoDz_Sub = (TH1F*)hRnPoDz_On->Clone();
	hRnPoDz_Sub->Sumw2();
	hRnPoDz_Sub->Add(hRnPoDz_Off,-1);
	hRnPoDz_Sub->Divide(hRnPoDz_Off);

	TH1F *hRate_Sub = (TH1F*)hRate_On->Clone();
	hRate_Sub->Sumw2();
	hRate_Sub->Add(hRate_Off,-1);
	hRate_Sub->Divide(hRate_Off);

	TH1F *hPoEnMean_Sub = (TH1F*)hPoEnMean_On->Clone();
	hPoEnMean_Sub->Sumw2();
	hPoEnMean_Sub->Add(hPoEnMean_Off,-1);

	TH1F *hPoEnSigma_Sub = (TH1F*)hPoEnSigma_On->Clone();
	hPoEnSigma_Sub->Sumw2();
	hPoEnSigma_Sub->Add(hPoEnSigma_Off,-1);

	TH1F *hRnPoDzMean_Sub = (TH1F*)hRnPoDzMean_On->Clone();
	hRnPoDzMean_Sub->Sumw2();
	hRnPoDzMean_Sub->Add(hRnPoDzMean_Off,-1);

	TH1F *hRnPoDzSigma_Sub = (TH1F*)hRnPoDzSigma_On->Clone();
	hRnPoDzSigma_Sub->Sumw2();
	hRnPoDzSigma_Sub->Add(hRnPoDzSigma_Off,-1);


	//-------------------------------------------------------------------------------------------------------
	TPaveText *pt;
	//int p_col = 8, d_col = 9;
	int p_col = 1, d_col = 2;
	int on_m = 21, off_m = 24;

cout<<"Plotting hists"<<endl;
	TLegend *leg;

	TCanvas *cRnPoEn = new TCanvas("cRnPoEn","RnPo En",700,900);
	cRnPoEn->Divide(1,2);
	cRnPoEn->cd(1);		
	gPad->SetGrid();
	hPoEn_On->GetYaxis()->SetRangeUser(0,3);
	hPoEn_On->GetYaxis()->SetTitle("Counts/5 keV/hr");
	hPoEn_On->SetMarkerColor(d_col);
	hPoEn_On->SetLineColor(d_col);
	hPoEn_On->SetMarkerStyle(on_m);
	hPoEn_On->Draw();
	hRnEn_On->SetMarkerColor(p_col);
	hRnEn_On->SetLineColor(p_col);
	hRnEn_On->SetMarkerStyle(on_m);
	hRnEn_On->Draw("same");
	hPoEn_Off->SetMarkerColor(d_col);
	hPoEn_Off->SetLineColor(d_col);
	hPoEn_Off->SetMarkerStyle(off_m);
	hPoEn_Off->Draw("same");
	hRnEn_Off->SetMarkerColor(p_col);
	hRnEn_Off->SetLineColor(p_col);
	hRnEn_Off->SetMarkerStyle(off_m);
	hRnEn_Off->Draw("same");
	leg = new TLegend(0.75,0.65,0.95,0.9);
	leg->AddEntry(hRnEn_On,"Rn^{219} Rx On","p");
	leg->AddEntry(hPoEn_On,"Po^{215} Rx On","p");
	leg->AddEntry(hRnEn_Off,"Rn^{219} Rx Off","p");
	leg->AddEntry(hPoEn_Off,"Po^{215} Rx Off","p");
	leg->Draw();
	cRnPoEn->cd(2);	
	gPad->SetGrid();
	hPoEn_Sub->GetYaxis()->SetRangeUser(-1,1);
	hPoEn_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hPoEn_Sub->SetMarkerColor(d_col);
	hPoEn_Sub->SetLineColor(d_col);
	hPoEn_Sub->SetMarkerStyle(2);	
	hPoEn_Sub->Draw();
	hRnEn_Sub->SetMarkerColor(p_col);
	hRnEn_Sub->SetLineColor(p_col);
	hRnEn_Sub->SetMarkerStyle(2);
	hRnEn_Sub->Draw("same");
	leg = new TLegend(0.83,0.7,0.95,0.9);
	leg->AddEntry(hRnEn_Sub,"Rn^{219}");
	leg->AddEntry(hPoEn_Sub,"Po^{215}");
	leg->Draw();
	cRnPoEn->SaveAs(Form("%s/RxOnVsRxOff/RnPoEn.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRnPoEn->SaveAs(Form("%s/RxOnVsRxOff/RnPoEn.png",gSystem->Getenv("AD_AC227_PLOTS")));
	

	TCanvas *cRnPoDz = new TCanvas("cRnPoDz","RnPo Dz",700,900);
	cRnPoDz->Divide(1,2);
	cRnPoDz->cd(1);
	gPad->SetGrid();
	hRnPoDz_On->GetYaxis()->SetRangeUser(0,1.5);
	hRnPoDz_On->GetYaxis()->SetTitle("Counts/0.25 cm/hr");
	hRnPoDz_On->SetMarkerColor(p_col);
	hRnPoDz_On->SetLineColor(p_col);
	hRnPoDz_On->SetMarkerStyle(on_m);
	hRnPoDz_On->Draw();
	hRnPoDz_Off->SetMarkerColor(d_col);
	hRnPoDz_Off->SetLineColor(d_col);
	hRnPoDz_Off->SetMarkerStyle(off_m);
	hRnPoDz_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hRnPoDz_On,"Rx On","p");
	leg->AddEntry(hRnPoDz_Off,"Rx Off","p");
	leg->Draw();		
	cRnPoDz->cd(2);
	gPad->SetGrid();
	hRnPoDz_Sub->GetYaxis()->SetRangeUser(-0.5,0.5);
	hRnPoDz_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hRnPoDz_Sub->SetMarkerColor(p_col);
	hRnPoDz_Sub->SetLineColor(p_col);
	hRnPoDz_Sub->SetMarkerStyle(2);
	hRnPoDz_Sub->Draw();
	cRnPoDz->SaveAs(Form("%s/RxOnVsRxOff/RnPoDz.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRnPoDz->SaveAs(Form("%s/RxOnVsRxOff/RnPoDz.png",gSystem->Getenv("AD_AC227_PLOTS")));
	
	TCanvas *cRate = new TCanvas("cRate","Rate Per Cell",700,900);
	cRate->Divide(1,2);
	cRate->cd(1);
	gPad->SetGrid();
	hRate_On->GetYaxis()->SetRangeUser(3,3.7);
	hRate_On->GetYaxis()->SetTitle("Rate [mHz]");
	hRate_On->GetXaxis()->SetTitle("Cell");
	hRate_On->SetMarkerColor(p_col);
	hRate_On->SetLineColor(p_col);
	hRate_On->SetMarkerStyle(on_m);
	hRate_On->Draw();
	hRate_Off->SetMarkerColor(d_col);
	hRate_Off->SetLineColor(d_col);
	hRate_Off->SetMarkerStyle(off_m);
	hRate_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hRate_On,"Rx On","p");
	leg->AddEntry(hRate_Off,"Rx Off","p");
	leg->Draw();		
	cRate->cd(2);
	gPad->SetGrid();
	hRate_Sub->GetYaxis()->SetRangeUser(-0.4,0.4);
	hRate_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hRate_Sub->GetXaxis()->SetTitle("Cell");
	hRate_Sub->SetMarkerColor(p_col);
	hRate_Sub->SetLineColor(p_col);
	hRate_Sub->SetMarkerStyle(2);
	hRate_Sub->Fit("pol0","0");
	hRate_Sub->Draw();
	cRate->SaveAs(Form("%s/RxOnVsRxOff/Rate.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RxOnVsRxOff/Rate.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cPoEnMean = new TCanvas("cPoEnMean","PoEnMean Per Cell",700,900);
	cPoEnMean->Divide(1,2);
	cPoEnMean->cd(1);
	gPad->SetGrid();
	hPoEnMean_On->GetYaxis()->SetRangeUser(793,813);
	hPoEnMean_On->GetYaxis()->SetTitle("E_{#mu Po} [keV]");
	hPoEnMean_On->GetXaxis()->SetTitle("Cell");
	hPoEnMean_On->SetMarkerColor(p_col);
	hPoEnMean_On->SetLineColor(p_col);
	hPoEnMean_On->SetMarkerStyle(on_m);
	hPoEnMean_On->Draw();
	hPoEnMean_Off->SetMarkerColor(d_col);
	hPoEnMean_Off->SetLineColor(d_col);
	hPoEnMean_Off->SetMarkerStyle(off_m);
	hPoEnMean_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hPoEnMean_On,"Rx On","p");
	leg->AddEntry(hPoEnMean_Off,"Rx Off","p");
	leg->Draw();		
	cPoEnMean->cd(2);
	gPad->SetGrid();
	hPoEnMean_Sub->GetYaxis()->SetRangeUser(-4,4);
	hPoEnMean_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hPoEnMean_Sub->GetXaxis()->SetTitle("Cell");
	hPoEnMean_Sub->SetMarkerColor(p_col);
	hPoEnMean_Sub->SetLineColor(p_col);
	hPoEnMean_Sub->SetMarkerStyle(2);
	hPoEnMean_Sub->Draw();
	hPoEnMean_Sub->Fit("pol0","0");
	cPoEnMean->SaveAs(Form("%s/RxOnVsRxOff/PoEnMean.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cPoEnMean->SaveAs(Form("%s/RxOnVsRxOff/PoEnMean.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cPoEnSigma = new TCanvas("cPoEnSigma","PoEnSigma Per Cell",700,900);
	cPoEnSigma->Divide(1,2);
	cPoEnSigma->cd(1);
	gPad->SetGrid();
	hPoEnSigma_On->GetYaxis()->SetRangeUser(37,53);
	hPoEnSigma_On->GetYaxis()->SetTitle("E_{#sigma Po} [keV]");
	hPoEnSigma_On->GetXaxis()->SetTitle("Cell");
	hPoEnSigma_On->SetMarkerColor(p_col);
	hPoEnSigma_On->SetLineColor(p_col);
	hPoEnSigma_On->SetMarkerStyle(on_m);
	hPoEnSigma_On->Draw();
	hPoEnSigma_Off->SetMarkerColor(d_col);
	hPoEnSigma_Off->SetLineColor(d_col);
	hPoEnSigma_Off->SetMarkerStyle(off_m);
	hPoEnSigma_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hPoEnSigma_On,"Rx On","p");
	leg->AddEntry(hPoEnSigma_Off,"Rx Off","p");
	leg->Draw();		
	cPoEnSigma->cd(2);
	gPad->SetGrid();
	hPoEnSigma_Sub->GetYaxis()->SetRangeUser(-4,4);
	hPoEnSigma_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hPoEnSigma_Sub->GetXaxis()->SetTitle("Cell");
	hPoEnSigma_Sub->SetMarkerColor(p_col);
	hPoEnSigma_Sub->SetLineColor(p_col);
	hPoEnSigma_Sub->SetMarkerStyle(2);
	hPoEnSigma_Sub->Fit("pol0","0");
	hPoEnSigma_Sub->Draw();
	cPoEnSigma->SaveAs(Form("%s/RxOnVsRxOff/PoEnSigma.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cPoEnSigma->SaveAs(Form("%s/RxOnVsRxOff/PoEnSigma.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRnPoDzMean = new TCanvas("cRnPoDzMean","RnPoDzMean Per Cell",700,900);
	cRnPoDzMean->Divide(1,2);
	cRnPoDzMean->cd(1);
	gPad->SetGrid();
	hRnPoDzMean_On->GetYaxis()->SetRangeUser(-3,3);
	hRnPoDzMean_On->GetYaxis()->SetTitle("dz_{#mu} [mm]");
	hRnPoDzMean_On->GetXaxis()->SetTitle("Cell");
	hRnPoDzMean_On->SetMarkerColor(p_col);
	hRnPoDzMean_On->SetLineColor(p_col);
	hRnPoDzMean_On->SetMarkerStyle(on_m);
	hRnPoDzMean_On->Draw();
	hRnPoDzMean_Off->SetMarkerColor(d_col);
	hRnPoDzMean_Off->SetLineColor(d_col);
	hRnPoDzMean_Off->SetMarkerStyle(off_m);
	hRnPoDzMean_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hRnPoDzMean_On,"Rx On","p");
	leg->AddEntry(hRnPoDzMean_Off,"Rx Off","p");
	leg->Draw();		
	cRnPoDzMean->cd(2);
	gPad->SetGrid();
	hRnPoDzMean_Sub->GetYaxis()->SetRangeUser(-5,5);
	hRnPoDzMean_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hRnPoDzMean_Sub->GetXaxis()->SetTitle("Cell");
	hRnPoDzMean_Sub->SetMarkerColor(p_col);
	hRnPoDzMean_Sub->SetLineColor(p_col);
	hRnPoDzMean_Sub->SetMarkerStyle(2);
	hRnPoDzMean_Sub->Fit("pol0","0");
	hRnPoDzMean_Sub->Draw();
	cRnPoDzMean->SaveAs(Form("%s/RxOnVsRxOff/RnPoDzMean.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRnPoDzMean->SaveAs(Form("%s/RxOnVsRxOff/RnPoDzMean.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRnPoDzSigma = new TCanvas("cRnPoDzSigma","RnPoDzSigma Per Cell",700,900);
	cRnPoDzSigma->Divide(1,2);
	cRnPoDzSigma->cd(1);
	gPad->SetGrid();
	hRnPoDzSigma_On->GetYaxis()->SetRangeUser(40,65);
	hRnPoDzSigma_On->GetYaxis()->SetTitle("dz_{#sigma} [mm]");
	hRnPoDzSigma_On->GetXaxis()->SetTitle("Cell");
	hRnPoDzSigma_On->SetMarkerColor(p_col);
	hRnPoDzSigma_On->SetLineColor(p_col);
	hRnPoDzSigma_On->SetMarkerStyle(on_m);
	hRnPoDzSigma_On->Draw();
	hRnPoDzSigma_Off->SetMarkerColor(d_col);
	hRnPoDzSigma_Off->SetLineColor(d_col);
	hRnPoDzSigma_Off->SetMarkerStyle(off_m);
	hRnPoDzSigma_Off->Draw("same");
	leg = new TLegend(0.8,0.75,0.95,0.9);
	leg->AddEntry(hRnPoDzSigma_On,"Rx On","p");
	leg->AddEntry(hRnPoDzSigma_Off,"Rx Off","p");
	leg->Draw();		
	cRnPoDzSigma->cd(2);
	gPad->SetGrid();
	hRnPoDzSigma_Sub->GetYaxis()->SetRangeUser(-10,10);
	hRnPoDzSigma_Sub->GetYaxis()->SetTitle("Rx On/Rx Off");
	hRnPoDzSigma_Sub->GetXaxis()->SetTitle("Cell");
	hRnPoDzSigma_Sub->SetMarkerColor(p_col);
	hRnPoDzSigma_Sub->SetLineColor(p_col);
	hRnPoDzSigma_Sub->SetMarkerStyle(2);
	hRnPoDzSigma_Sub->Fit("pol0","0");
	hRnPoDzSigma_Sub->Draw();
	cRnPoDzSigma->SaveAs(Form("%s/RxOnVsRxOff/RnPoDzSigma.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRnPoDzSigma->SaveAs(Form("%s/RxOnVsRxOff/RnPoDzSigma.png",gSystem->Getenv("AD_AC227_PLOTS")));


	return 0;
} 	//end PlotDistributionsVsTime

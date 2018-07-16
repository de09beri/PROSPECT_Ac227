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

const int NUMCOLS = 14;

TGraphErrors *makeRelGr(TGraphErrors *gr, double mean, double meanErr){
	TGraphErrors *grRel = (TGraphErrors*)gr->Clone();
	int numPt = gr->GetN();

	double rel, relErr; 
	double grx, gry, gryErr;

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		rel = gry/mean;
		relErr = rel * sqrt(pow(gryErr/gry,2) + pow(meanErr/mean,2));	

		grRel->SetPoint(i,grx,rel);
		grRel->SetPointError(i,0,relErr);
	}

	return grRel;
}


void PlotRnPoVsTime(){

	setup_PROSPECT_style();
    gROOT->ForceStyle();

	TGraphErrors *grRate[NUMCOLS],	 	   *grRelRate[NUMCOLS];
	TGraphErrors *grPoPSDMean[NUMCOLS],    *grRelPoPSDMean[NUMCOLS];
	TGraphErrors *grPoPSDSigma[NUMCOLS],   *grRelPoPSDSigma[NUMCOLS];
	TGraphErrors *grPoEnMean[NUMCOLS],     *grRelPoEnMean[NUMCOLS];
	TGraphErrors *grPoEnSigma[NUMCOLS],    *grRelPoEnSigma[NUMCOLS];
	TGraphErrors *grPoPosSigma[NUMCOLS],   *grRelPoPosSigma[NUMCOLS];
	TGraphErrors *grRnPoDzSigma[NUMCOLS],  *grRelRnPoDzSigma[NUMCOLS];
	

	TFile *f = new TFile(Form("%s/Ac227_ColGraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	for(int i=0;i<NUMCOLS;i++){
		grRate[i]        = (TGraphErrors*)f->Get(Form("grRate_%i",i));
		grPoPSDMean[i]   = (TGraphErrors*)f->Get(Form("grPoPSDMean_%i",i));
		grPoPSDSigma[i]  = (TGraphErrors*)f->Get(Form("grPoPSDSigma_%i",i));
		grPoEnMean[i]    = (TGraphErrors*)f->Get(Form("grPoEnMean_%i",i));
		grPoEnSigma[i]   = (TGraphErrors*)f->Get(Form("grPoEnSigma_%i",i));
		grPoPosSigma[i]  = (TGraphErrors*)f->Get(Form("grPoPosSigma_%i",i));
		grRnPoDzSigma[i] = (TGraphErrors*)f->Get(Form("grRnPoDzSigma_%i",i));
	}

	f->Close();
	
	//-------------------------------------------------------------------------------------------------------
	for(int i=0;i<NUMCOLS;i++){
		grRate[i]->Fit("pol0","Q0");
		grRelRate[i] = makeRelGr(grRate[i],grRate[i]->GetFunction("pol0")->GetParameter(0),grRate[i]->GetFunction("pol0")->GetParError(0));

		grPoPSDMean[i]->Fit("pol0","Q0");
		grRelPoPSDMean[i] = makeRelGr(grPoPSDMean[i], grPoPSDMean[i]->GetFunction("pol0")->GetParameter(0), grPoPSDMean[i]->GetFunction("pol0")->GetParError(0));

		grPoPSDSigma[i]->Fit("pol0","Q0");
		grRelPoPSDSigma[i] = makeRelGr(grPoPSDSigma[i], grPoPSDSigma[i]->GetFunction("pol0")->GetParameter(0), grPoPSDSigma[i]->GetFunction("pol0")->GetParError(0));

		grPoEnMean[i]->Fit("pol0","Q0");
		grRelPoEnMean[i] = makeRelGr(grPoEnMean[i], grPoEnMean[i]->GetFunction("pol0")->GetParameter(0), grPoEnMean[i]->GetFunction("pol0")->GetParError(0));

		grPoEnSigma[i]->Fit("pol0","Q0");
		grRelPoEnSigma[i] = makeRelGr(grPoEnSigma[i], grPoEnSigma[i]->GetFunction("pol0")->GetParameter(0), grPoEnSigma[i]->GetFunction("pol0")->GetParError(0));

		grPoPosSigma[i]->Fit("pol0","Q0");
		grRelPoPosSigma[i] = makeRelGr(grPoPosSigma[i], grPoPosSigma[i]->GetFunction("pol0")->GetParameter(0), grPoPosSigma[i]->GetFunction("pol0")->GetParError(0));

		grRnPoDzSigma[i]->Fit("pol0","Q0");
		grRelRnPoDzSigma[i] = makeRelGr(grRnPoDzSigma[i], grRnPoDzSigma[i]->GetFunction("pol0")->GetParameter(0), grRnPoDzSigma[i]->GetFunction("pol0")->GetParError(0));
	}

	//-------------------------------------------------------------------------------------------------------
	int col[NUMCOLS] = {46,41,40,36,30,28,11,9,8,7,6,4,2,1};
	const char *xLabel = "Time";
	TLegend *leg;

	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate vs time",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelRate[0]->SetLineColor(col[0]);
	grRelRate[0]->SetMarkerColor(col[0]);
	grRelRate[0]->GetXaxis()->SetTitle(xLabel);
	grRelRate[0]->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");
	grRelRate[0]->GetXaxis()->SetTimeDisplay(1);
	grRelRate[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRate[0]->GetYaxis()->SetRangeUser(0.75,1.15);
	grRelRate[0]->Draw("ALP");
	leg->AddEntry(grRelRate[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelRate[i]->SetLineColor(col[i]);
	grRelRate[i]->SetMarkerColor(col[i]);
	grRelRate[i]->Draw("LP");
	leg->AddEntry(grRelRate[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelRate->SaveAs(Form("%s/ColRelativeRateVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/ColRelativeRateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelPoPSDMean[0]->SetLineColor(col[0]);
	grRelPoPSDMean[0]->SetMarkerColor(col[0]);
	grRelPoPSDMean[0]->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean[0]->GetYaxis()->SetTitle("PSD_{#muPo}/#LTPSD_{#muPo}#GT");
	grRelPoPSDMean[0]->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDMean[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDMean[0]->GetYaxis()->SetRangeUser(0.97,1.04);
	grRelPoPSDMean[0]->Draw("ALP");
	leg->AddEntry(grRelPoPSDMean[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelPoPSDMean[i]->SetLineColor(col[i]);
	grRelPoPSDMean[i]->SetMarkerColor(col[i]);
	grRelPoPSDMean[i]->Draw("LP");
	leg->AddEntry(grRelPoPSDMean[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelPoPSDMean->SaveAs(Form("%s/ColRelativePoPSDMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/ColRelativePoPSDMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelPoPSDSigma[0]->SetLineColor(col[0]);
	grRelPoPSDSigma[0]->SetMarkerColor(col[0]);
	grRelPoPSDSigma[0]->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma[0]->GetYaxis()->SetTitle("PSD_{#sigmaPo}/#LTPSD_{#sigmaPo}#GT");
	grRelPoPSDSigma[0]->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDSigma[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDSigma[0]->GetYaxis()->SetRangeUser(0.92,1.1);
	grRelPoPSDSigma[0]->Draw("ALP");
	leg->AddEntry(grRelPoPSDSigma[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelPoPSDSigma[i]->SetLineColor(col[i]);
	grRelPoPSDSigma[i]->SetMarkerColor(col[i]);
	grRelPoPSDSigma[i]->Draw("LP");
	leg->AddEntry(grRelPoPSDSigma[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelPoPSDSigma->SaveAs(Form("%s/ColRelativePoPSDSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/ColRelativePoPSDSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po En Mean",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelPoEnMean[0]->SetLineColor(col[0]);
	grRelPoEnMean[0]->SetMarkerColor(col[0]);
	grRelPoEnMean[0]->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean[0]->GetYaxis()->SetTitle("E_{#muPo}/#LTE_{#muPo}#GT");
	grRelPoEnMean[0]->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnMean[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnMean[0]->GetYaxis()->SetRangeUser(0.99,1.01);
	grRelPoEnMean[0]->Draw("ALP");
	leg->AddEntry(grRelPoEnMean[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelPoEnMean[i]->SetLineColor(col[i]);
	grRelPoEnMean[i]->SetMarkerColor(col[i]);
	grRelPoEnMean[i]->Draw("LP");
	leg->AddEntry(grRelPoEnMean[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelPoEnMean->SaveAs(Form("%s/ColRelativePoEnMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/ColRelativePoEnMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po En Sigma",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelPoEnSigma[0]->SetLineColor(col[0]);
	grRelPoEnSigma[0]->SetMarkerColor(col[0]);
	grRelPoEnSigma[0]->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma[0]->GetYaxis()->SetTitle("E_{#sigmaPo}/#LTE_{#sigmaPo}#GT");
	grRelPoEnSigma[0]->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnSigma[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnSigma[0]->GetYaxis()->SetRangeUser(0.8,1.15);
	grRelPoEnSigma[0]->Draw("ALP");
	leg->AddEntry(grRelPoEnSigma[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelPoEnSigma[i]->SetLineColor(col[i]);
	grRelPoEnSigma[i]->SetMarkerColor(col[i]);
	grRelPoEnSigma[i]->Draw("LP");
	leg->AddEntry(grRelPoEnSigma[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelPoEnSigma->SaveAs(Form("%s/ColRelativePoEnSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/ColRelativePoEnSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Pos Sigma",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelPoPosSigma[0]->SetLineColor(col[0]);
	grRelPoPosSigma[0]->SetMarkerColor(col[0]);
	grRelPoPosSigma[0]->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma[0]->GetYaxis()->SetTitle("z_{#sigmaPo}/#LTz_{#sigmaPo}#GT");
	grRelPoPosSigma[0]->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosSigma[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosSigma[0]->GetYaxis()->SetRangeUser(0.85,1.15);
	grRelPoPosSigma[0]->Draw("ALP");
	leg->AddEntry(grRelPoPosSigma[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelPoPosSigma[i]->SetLineColor(col[i]);
	grRelPoPosSigma[i]->SetMarkerColor(col[i]);
	grRelPoPosSigma[i]->Draw("LP");
	leg->AddEntry(grRelPoPosSigma[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelPoPosSigma->SaveAs(Form("%s/ColRelativePoPosSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/ColRelativePoPosSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	leg = new TLegend(0.94,0.30,0.999,0.98);
	grRelRnPoDzSigma[0]->SetLineColor(col[0]);
	grRelRnPoDzSigma[0]->SetMarkerColor(col[0]);
	grRelRnPoDzSigma[0]->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma[0]->GetYaxis()->SetTitle("dz_{#sigmaRnPo}/#LTdz_{#sigmaRnPo}#GT");
	grRelRnPoDzSigma[0]->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzSigma[0]->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzSigma[0]->GetYaxis()->SetRangeUser(0.8,1.2);
	grRelRnPoDzSigma[0]->Draw("ALP");
	leg->AddEntry(grRelRnPoDzSigma[0],"0","l");
	for(int i=1;i<NUMCOLS;i++){
	grRelRnPoDzSigma[i]->SetLineColor(col[i]);
	grRelRnPoDzSigma[i]->SetMarkerColor(col[i]);
	grRelRnPoDzSigma[i]->Draw("LP");
	leg->AddEntry(grRelRnPoDzSigma[i],Form("%i",i),"l");
	}
	leg->Draw();
	cRelRnPoDzSigma->SaveAs(Form("%s/ColRelativeRnPoDzSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/ColRelativeRnPoDzSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

} 	//end PlotRnPoVsTime

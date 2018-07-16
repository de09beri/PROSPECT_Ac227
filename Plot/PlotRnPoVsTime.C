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

	TFile *f = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")));

	TGraphErrors *grRate 		= (TGraphErrors*)f->Get("grRate");
	TGraphErrors *grTotEff 		= (TGraphErrors*)f->Get("grTotEff");
	TGraphErrors *grLifetime 	= (TGraphErrors*)f->Get("grLifetime");
	TGraphErrors *grPoPSDMean 	= (TGraphErrors*)f->Get("grPoPSDMean");
	TGraphErrors *grPoPSDSigma 	= (TGraphErrors*)f->Get("grPoPSDSigma");
	TGraphErrors *grPoEnMean 	= (TGraphErrors*)f->Get("grPoEnMean");
	TGraphErrors *grPoEnSigma 	= (TGraphErrors*)f->Get("grPoEnSigma");
	TGraphErrors *grPoPosMean 	= (TGraphErrors*)f->Get("grPoPosMean");
	TGraphErrors *grPoPosSigma 	= (TGraphErrors*)f->Get("grPoPosSigma");
	TGraphErrors *grRnPoDzMean 	= (TGraphErrors*)f->Get("grRnPoDzMean");
	TGraphErrors *grRnPoDzSigma = (TGraphErrors*)f->Get("grRnPoDzSigma");

	f->Close();
	
	//-------------------------------------------------------------------------------------------------------
//	grRate->RemovePoint(10);
//	grRate->RemovePoint(1);
cout<<"\n RATE"<<endl;
	grRate->Fit("pol0","0");
	TGraphErrors *grRelRate = makeRelGr(grRate, grRate->GetFunction("pol0")->GetParameter(0), grRate->GetFunction("pol0")->GetParError(0));

cout<<"\n PSD MEAN"<<endl;
	grPoPSDMean->Fit("pol0","0");
	TGraphErrors *grRelPoPSDMean = makeRelGr(grPoPSDMean, grPoPSDMean->GetFunction("pol0")->GetParameter(0), grPoPSDMean->GetFunction("pol0")->GetParError(0));

cout<<"\n PSD SIGMA"<<endl;
	grPoPSDSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPSDSigma = makeRelGr(grPoPSDSigma, grPoPSDSigma->GetFunction("pol0")->GetParameter(0), grPoPSDSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n ENERGY MEAN"<<endl;
	grPoEnMean->Fit("pol0","0");
	TGraphErrors *grRelPoEnMean = makeRelGr(grPoEnMean, grPoEnMean->GetFunction("pol0")->GetParameter(0), grPoEnMean->GetFunction("pol0")->GetParError(0));

cout<<"\n ENERGY SIGMA"<<endl;
	grPoEnSigma->Fit("pol0","0");
	TGraphErrors *grRelPoEnSigma = makeRelGr(grPoEnSigma, grPoEnSigma->GetFunction("pol0")->GetParameter(0), grPoEnSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n POSITION MEAN"<<endl;
	grPoPosMean->Fit("pol0","0");
	TGraphErrors *grRelPoPosMean = makeRelGr(grPoPosMean, grPoPosMean->GetFunction("pol0")->GetParameter(0), grPoPosMean->GetFunction("pol0")->GetParError(0));

cout<<"\n POSITION RMS"<<endl;
	grPoPosSigma->Fit("pol0","0");
	TGraphErrors *grRelPoPosSigma = makeRelGr(grPoPosSigma, grPoPosSigma->GetFunction("pol0")->GetParameter(0), grPoPosSigma->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ MEAN"<<endl;
	grRnPoDzMean->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzMean = makeRelGr(grRnPoDzMean, grRnPoDzMean->GetFunction("pol0")->GetParameter(0), grRnPoDzMean->GetFunction("pol0")->GetParError(0));

cout<<"\n DZ SIGMA"<<endl;
	grRnPoDzSigma->Fit("pol0","0");
	TGraphErrors *grRelRnPoDzSigma = makeRelGr(grRnPoDzSigma, grRnPoDzSigma->GetFunction("pol0")->GetParameter(0), grRnPoDzSigma->GetFunction("pol0")->GetParError(0));
	
	//-------------------------------------------------------------------------------------------------------
	const char *xLabel = "Time";

	TCanvas *cRate = new TCanvas("cRate","Rate vs time",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	cRate->SaveAs(Form("%s/RateVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRate->SaveAs(Form("%s/RateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRate = new TCanvas("cRelRate","Relative rate vs time",1000,400);
	grRelRate->GetXaxis()->SetTitle(xLabel);
	grRelRate->GetYaxis()->SetTitle("R_{RnPo}/#LTR_{RnPo}#GT");  
	grRelRate->GetXaxis()->SetTimeDisplay(1);
	grRelRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRate->Draw("AP");
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRate->SaveAs(Form("%s/RelativeRateVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


	int numPt = grRate->GetN();
	double grxStart, grxEnd, gry;
	grRate->GetPoint(0,grxStart,gry);
	grRate->GetPoint(numPt-1,grxEnd,gry);

	//TF1 *fAcExp = new TF1("fAcExp","[0]*exp(-((x-1514782800)*log(2)/[1]))",grxStart,grxEnd);
	TF1 *fAcExp = new TF1("fAcExp","[0]*exp(-(x*log(2)/[1]))",grxStart,grxEnd);
	fAcExp->SetParameters(0.40,684300000);

	TCanvas *cRateFit = new TCanvas("cRate","Rate vs time Fit",1000,400);
	grRate->GetXaxis()->SetTitle(xLabel);
	grRate->GetYaxis()->SetTitle("R_{RnPo} [Hz]");
	grRate->GetXaxis()->SetTimeDisplay(1);
	grRate->GetXaxis()->SetTimeFormat("%m/%d");
	grRate->Draw("AP");
	grRate->Fit(fAcExp,"R0");
	fAcExp->Draw("same");
	TPaveText *pv = new TPaveText(0.4,0.8,0.6,0.98,"NDC");
	pv->AddText(Form("#Chi^{2}/NDF  %f/%d",fAcExp->GetChisquare(),fAcExp->GetNDF()));
	pv->AddText(Form("Prob  %f",fAcExp->GetProb()));
	pv->AddText(Form("R_{0}   %.2f #pm %.2f Hz",fAcExp->GetParameter(0),fAcExp->GetParError(0)));
	pv->AddText(Form("t_{1/2}   %.2f #pm %.2f yrs",fAcExp->GetParameter(1)/(60.0*60.0*24.0*365.0),fAcExp->GetParError(1)/(60.0*60.0*24.0*365.0)));
	pv->Draw();
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRateFit->SaveAs(Form("%s/RateVsTime_Fit.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDMean = new TCanvas("cRelPoPSDMean","Relative Po PSD Mean",1000,400);
	grRelPoPSDMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDMean->GetYaxis()->SetTitle("PSD_{#muPo}/#LTPSD_{#muPo}#GT");
	grRelPoPSDMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDMean->Draw("AP");
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDMean->SaveAs(Form("%s/RelativePoPSDMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPSDSigma = new TCanvas("cRelPoPSDSigma","Relative Po PSD Sigma",1000,400);
	grRelPoPSDSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPSDSigma->GetYaxis()->SetTitle("PSD_{#sigmaPo}/#LTPSD_{#sigmaPo}#GT");
	grRelPoPSDSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPSDSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPSDSigma->Draw("AP");
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPSDSigma->SaveAs(Form("%s/RelativePoPSDSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnMean = new TCanvas("cRelPoEnMean","Relative Po Energy Mean",1000,400);
	grRelPoEnMean->GetXaxis()->SetTitle(xLabel);
	grRelPoEnMean->GetYaxis()->SetTitle("E_{#muPo}/#LTE_{#muPo}#GT");
	grRelPoEnMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnMean->Draw("AP");
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnMean->SaveAs(Form("%s/RelativePoEnMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoEnSigma = new TCanvas("cRelPoEnSigma","Relative Po Energy Sigma",1000,400);
	grRelPoEnSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoEnSigma->GetYaxis()->SetTitle("E_{#sigmaPo}/#LTE_{#sigmaPo}#GT");
	grRelPoEnSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoEnSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoEnSigma->Draw("AP");
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoEnSigma->SaveAs(Form("%s/RelativePoEnSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosMean = new TCanvas("cRelPoPosMean","Relative Po Position Mean",1000,400);
	grRelPoPosMean->GetXaxis()->SetTitle(xLabel);
	grRelPoPosMean->GetYaxis()->SetTitle("z_{#muPo}/#LTz_{#muPo}#GT");
	grRelPoPosMean->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosMean->Draw("AP");
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosMean->SaveAs(Form("%s/RelativePoPosMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelPoPosSigma = new TCanvas("cRelPoPosSigma","Relative Po Position Sigma",1000,400);
	grRelPoPosSigma->GetXaxis()->SetTitle(xLabel);
	grRelPoPosSigma->GetYaxis()->SetTitle("z_{RMS Po}/#LTz_{RMS Po}#GT");
	grRelPoPosSigma->GetXaxis()->SetTimeDisplay(1);
	grRelPoPosSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelPoPosSigma->Draw("AP");
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelPoPosSigma->SaveAs(Form("%s/RelativePoPosSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzMean = new TCanvas("cRelRnPoDzMean","Relative RnPo Dz Mean",1000,400);
	grRelRnPoDzMean->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzMean->GetYaxis()->SetTitle("dz_{#mu}/#LTdz_{#mu}#GT");
	grRelRnPoDzMean->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzMean->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzMean->Draw("AP");
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzMean->SaveAs(Form("%s/RelativeRnPoDzMeanVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));

	TCanvas *cRelRnPoDzSigma = new TCanvas("cRelRnPoDzSigma","Relative RnPo Dz Sigma",1000,400);
	grRelRnPoDzSigma->GetXaxis()->SetTitle(xLabel);
	grRelRnPoDzSigma->GetYaxis()->SetTitle("dz_{#sigma}/#LTdz_{#sigma}#GT");
	grRelRnPoDzSigma->GetXaxis()->SetTimeDisplay(1);
	grRelRnPoDzSigma->GetXaxis()->SetTimeFormat("%m/%d");
	grRelRnPoDzSigma->Draw("AP");
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
	cRelRnPoDzSigma->SaveAs(Form("%s/RelativeRnPoDzSigmaVsTime.png",gSystem->Getenv("AD_AC227_PLOTS")));


} 	//end PlotRnPoVsTime

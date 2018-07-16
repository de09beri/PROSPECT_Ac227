#include "RNPO.C"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TMath.h"
#include "TLatex.h"
#include "TVectorD.h"
#include "TLegend.h"

using namespace std;

void PlotCell136(){

double TIMEOFFSET = 10.0*2.57;

RNPO *rnpo = new RNPO();

TH1F *hSelectRnPSD_Bump = new TH1F("hSelectRnPSD_Bump","",100,0.18,0.36);
TH1F *hBGRnPSD_Bump = new TH1F("hBGRnPSD_Bump","",100,0.18,0.36);

TH1F *hSelectRnPSD_noBump = new TH1F("hSelectRnPSD_noBump","",100,0.18,0.36);
TH1F *hBGRnPSD_noBump = new TH1F("hBGRnPSD_noBump","",100,0.18,0.36);

TH1F *hSelectPoPSD_Bump = new TH1F("hSelectPoPSD_Bump","",100,0.18,0.36);
TH1F *hBGPoPSD_Bump = new TH1F("hBGPoPSD_Bump","",100,0.18,0.36);

TH1F *hSelectPoPSD_noBump = new TH1F("hSelectPoPSD_noBump","",100,0.18,0.36);
TH1F *hBGPoPSD_noBump = new TH1F("hBGPoPSD_noBump","",100,0.18,0.36);

TH1F *hSelectDt_Bump = new TH1F("hSelectDt_Bump","",100,0,12.85);
TH1F *hBGDt_Bump = new TH1F("hBGDt_Bump","",100,0,12.85);

TH1F *hSelectDt_noBump = new TH1F("hSelectDt_noBump","",100,0,12.85);
TH1F *hBGDt_noBump = new TH1F("hBGDt_noBump","",100,0,12.85);


rnpo->fChain->Draw("p_PSD>>hSelectRnPSD_Bump","p_seg==136 && p_z>-250 && p_z<100");
rnpo->fChain->Draw("f_PSD>>hBGRnPSD_Bump","f_seg==136 && f_z>-250 && f_z<100");

rnpo->fChain->Draw("p_PSD>>hSelectRnPSD_noBump","p_seg==136 && p_z>100 && p_z<450");
rnpo->fChain->Draw("f_PSD>>hBGRnPSD_noBump","f_seg==136 && f_z>100 && f_z<450");

rnpo->fChain->Draw("d_PSD>>hSelectPoPSD_Bump","p_seg==136 && d_z>-250 && d_z<100");
rnpo->fChain->Draw("d_PSD>>hBGPoPSD_Bump","f_seg==136 && d_z>-250 && d_z<100");

rnpo->fChain->Draw("d_PSD>>hSelectPoPSD_noBump","p_seg==136 && d_z>100 && d_z<450");
rnpo->fChain->Draw("d_PSD>>hBGPoPSD_noBump","f_seg==136 && d_z>100 && d_z<450");


rnpo->fChain->Draw("(d_t - p_t)*(1e-6)>>hSelectDt_Bump","d_seg==136 && ((d_z>-250 && d_z<100) || (p_z>-250 && p_z<100) || (f_z>-250 && f_z<100))"); 
rnpo->fChain->Draw("((f_t - d_t)*(1e-6)-(10.0*2.57))>>hBGDt_Bump","d_seg==136 && ((d_z>-250 && d_z<100) || (p_z>-250 && p_z<100) || (f_z>-250 && f_z<100))"); 

rnpo->fChain->Draw("(d_t - p_t)*(1e-6)>>hSelectDt_noBump","d_seg==136 && ((d_z>100 && d_z<450) || (p_z>100 && p_z<450) || (f_z>100 && f_z<450))"); 
rnpo->fChain->Draw("((f_t - d_t)*(1e-6)-(10.0*2.57))>>hBGDt_noBump","d_seg==136 && ((d_z>100 && d_z<450) || (p_z>100 && p_z<450) || (f_z>100 && f_z<450))"); 

//=============================================================

TH1F *hRnPSD_Bump = (TH1F*)hSelectRnPSD_Bump->Clone();
hRnPSD_Bump->SetName("hRnPSD_Bump");
hRnPSD_Bump->Sumw2();
hRnPSD_Bump->Add(hBGRnPSD_Bump,-1);

TH1F *hRnPSD_noBump = (TH1F*)hSelectRnPSD_noBump->Clone();
hRnPSD_noBump->SetName("hRnPSD_noBump");
hRnPSD_noBump->Sumw2();
hRnPSD_noBump->Add(hBGRnPSD_noBump,-1);

TH1F *hPoPSD_Bump = (TH1F*)hSelectPoPSD_Bump->Clone();
hPoPSD_Bump->SetName("hPoPSD_Bump");
hPoPSD_Bump->Sumw2();
hPoPSD_Bump->Add(hBGPoPSD_Bump,-1);

TH1F *hPoPSD_noBump = (TH1F*)hSelectPoPSD_noBump->Clone();
hPoPSD_noBump->SetName("hPoPSD_noBump");
hPoPSD_noBump->Sumw2();
hPoPSD_noBump->Add(hBGPoPSD_noBump,-1);

TH1F *hRnPSD_Sub = (TH1F*)hRnPSD_Bump->Clone();
hRnPSD_Sub->SetName("hRnPSD_Sub");
hRnPSD_Sub->Sumw2();
hRnPSD_Sub->Add(hRnPSD_noBump,-1);

TH1F *hPoPSD_Sub = (TH1F*)hPoPSD_Bump->Clone();
hPoPSD_Sub->SetName("hPoPSD_Sub");
hPoPSD_Sub->Sumw2();
hPoPSD_Sub->Add(hPoPSD_noBump,-1);


TH1F *hRnPoDt_Bump = (TH1F*)hSelectDt_Bump->Clone();
hRnPoDt_Bump->SetName("hRnPoDt_Bump");
hRnPoDt_Bump->Sumw2();
hRnPoDt_Bump->Add(hBGDt_Bump,-1);

TH1F *hRnPoDt_noBump = (TH1F*)hSelectDt_noBump->Clone();
hRnPoDt_noBump->SetName("hRnPoDt_noBump");
hRnPoDt_noBump->Sumw2();
hRnPoDt_noBump->Add(hBGDt_noBump,-1);

TH1F *hRnPoDt_Sub = (TH1F*)hRnPoDt_Bump->Clone();
hRnPoDt_Sub->SetName("hRnPoDt_Sub");
hRnPoDt_Sub->Sumw2();
hRnPoDt_Sub->Add(hRnPoDt_noBump,-1);

gStyle->SetOptFit(1111);
TLegend *leg;


TCanvas *cPSD = new TCanvas("cPSD","PSD",1400,900);
cPSD->Divide(2,2);
cPSD->cd(1);
gPad->SetGrid();
hRnPSD_Bump->SetTitle(";Rn^{219} PSD;Counts");
hRnPSD_Bump->SetLineWidth(2);
hRnPSD_Bump->Draw();
hRnPSD_noBump->SetLineWidth(2);
hRnPSD_noBump->SetLineColor(kRed);
hRnPSD_noBump->SetMarkerColor(kRed);
hRnPSD_noBump->Draw("same");
leg = new TLegend(0.7,0.8,0.98,0.98);
leg->AddEntry(hRnPSD_Bump,"-250<z<100 mm","l");
leg->AddEntry(hRnPSD_noBump,"100<z<450 mm","l");
leg->Draw();
cPSD->cd(2);
gPad->SetGrid();
hPoPSD_Bump->SetTitle(";Po^{215} PSD;Counts");
hPoPSD_Bump->SetLineWidth(2);
hPoPSD_Bump->Draw();
hPoPSD_noBump->SetLineWidth(2);
hPoPSD_noBump->SetLineColor(kRed);
hPoPSD_noBump->SetMarkerColor(kRed);
hPoPSD_noBump->Draw("same");
leg = new TLegend(0.7,0.8,0.98,0.98);
leg->AddEntry(hPoPSD_Bump,"-250<z<100 mm","l");
leg->AddEntry(hPoPSD_noBump,"100<z<450 mm","l");
leg->Draw();
cPSD->cd(3);
gPad->SetGrid();
hRnPSD_Sub->SetTitle(";PSD_z(-250,100) - PSD_z(100,450);Counts"); 
hRnPSD_Sub->SetLineWidth(2);
hRnPSD_Sub->SetLineColor(kBlack);
hRnPSD_Sub->SetMarkerColor(kBlack);
hRnPSD_Sub->GetYaxis()->SetRangeUser(-30,80);
hRnPSD_Sub->Draw();
hRnPSD_Sub->Fit("pol0");
cPSD->cd(4);
gPad->SetGrid();
hPoPSD_Sub->SetTitle(";PSD_z(-250,100) - PSD_z(100,450);Counts"); 
hPoPSD_Sub->SetLineWidth(2);
hPoPSD_Sub->SetLineColor(kBlack);
hPoPSD_Sub->SetMarkerColor(kBlack);
hPoPSD_Sub->GetYaxis()->SetRangeUser(-30,80);
hPoPSD_Sub->Draw();
hPoPSD_Sub->Fit("pol0");
cPSD->SaveAs(Form("%s/SubtractedPSD_Cell136.png",gSystem->Getenv("AD_AC227_PLOTS")));


TCanvas *cDt = new TCanvas("cDt","Dt",700,900);
cDt->Divide(1,2);
cDt->cd(1);
gPad->SetGrid();
hRnPoDt_Bump->SetTitle("; #Delta t [ms];Counts");
hRnPoDt_Bump->SetLineWidth(2);
hRnPoDt_Bump->Draw();
hRnPoDt_noBump->SetLineWidth(2);
hRnPoDt_noBump->SetLineColor(kRed);
hRnPoDt_noBump->SetMarkerColor(kRed);
hRnPoDt_noBump->Draw("same");
leg = new TLegend(0.7,0.8,0.98,0.98);
leg->AddEntry(hRnPoDt_Bump,"-250<z<100 mm","l");
leg->AddEntry(hRnPoDt_noBump,"100<z<450 mm","l");
leg->Draw();
cDt->cd(2);
gPad->SetGrid();
hRnPoDt_Sub->SetTitle(";#Delta t_z(-250,100) - #Delta t_z(100,450);Counts"); 
hRnPoDt_Sub->SetLineWidth(2);
hRnPoDt_Sub->SetLineColor(kBlack);
hRnPoDt_Sub->SetMarkerColor(kBlack);
hRnPoDt_Sub->Draw();
TF1 *fExp = new TF1("fExp","[0]*exp(-(x*log(2))/[1]) + [2]",0,12.85);
fExp->SetParameters(0.003,4.3,0.002);
hRnPoDt_Sub->Fit("fExp","R");
cDt->SaveAs(Form("%s/SubtractedDt_Cell136.png",gSystem->Getenv("AD_AC227_PLOTS")));


}

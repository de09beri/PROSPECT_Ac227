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

using namespace std;

void PlotTime(){

double TIMEOFFSET = 10.0*2.57;

RNPO *rnpo = new RNPO();

TH1F *hSelectDt = new TH1F("hSelectDt","",150,0,0.5);
TH1F *hBGDt = new TH1F("hBGDt","",150,0,0.5);

rnpo->fChain->Draw("(d_t - p_t)*(1e-6)>>hSelectDt","p_seg>-1 && d_E>0.66"); 
rnpo->fChain->Draw("((f_t - d_t)*(1e-6)-(10.0*2.57))>>hBGDt","f_seg>-1 && d_E>0.66"); 

TH1F *hRnPoDt = (TH1F*)hSelectDt->Clone();
hRnPoDt->SetName("hRnPoDt");
hRnPoDt->Sumw2();
hRnPoDt->Add(hBGDt,-1);

double binwidth = 0.5/150.0;
TF1 *f = new TF1("f","[0]*exp(-log(2)*x/[1]) + [2]*exp(-log(2)*x/[3]) + [4]",binwidth*2,0.5);
f->SetParameters(3000,1.78,5000,0.005,1000);

gStyle->SetOptFit(1111);

TCanvas *c = new TCanvas("c","c",1);
hRnPoDt->Draw();
hRnPoDt->Fit(f,"R");

TFile *fout = new TFile("DtHist.root","RECREATE");
hRnPoDt->Write();
fout->Close();

}

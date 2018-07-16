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

void PlotEnergy(){

RNPO *rnpo = new RNPO();

TH1F *h_pE_mult1 = new TH1F("h_pE_mult1","",150,0.5,1.15);
TH1F *h_fE_mult1 = new TH1F("h_fE_mult1","",150,0.5,1.15);

TH1F *h_pE_mult2 = new TH1F("h_pE_mult2","",150,0.5,1.15);
TH1F *h_fE_mult2 = new TH1F("h_fE_mult2","",150,0.5,1.15);

TH1F *h_pE_mult3 = new TH1F("h_pE_mult3","",150,0.5,1.15);
TH1F *h_fE_mult3 = new TH1F("h_fE_mult3","",150,0.5,1.15);

TH1F *h_pE_mult4 = new TH1F("h_pE_mult4","",150,0.5,1.15);
TH1F *h_fE_mult4 = new TH1F("h_fE_mult4","",150,0.5,1.15);

TH1F *h_pE_mult5 = new TH1F("h_pE_mult5","",150,0.5,1.15);
TH1F *h_fE_mult5 = new TH1F("h_fE_mult5","",150,0.5,1.15);


rnpo->fChain->Draw("p_Etot>>h_pE_mult1","p_seg>0 && p_clustMult==1");
rnpo->fChain->Draw("f_Etot>>h_fE_mult1","f_seg>0 && f_clustMult==1","same");
h_pE_mult1->Sumw2();
TH1F *h_RnE_mult1 = (TH1F*)h_pE_mult1->Clone();
h_RnE_mult1->Add(h_fE_mult1,-1);

rnpo->fChain->Draw("p_Etot>>h_pE_mult2","p_seg>0 && p_clustMult==2");
rnpo->fChain->Draw("f_Etot>>h_fE_mult2","f_seg>0 && f_clustMult==2","same");
h_pE_mult2->Sumw2();
TH1F *h_RnE_mult2 = (TH1F*)h_pE_mult2->Clone();
h_RnE_mult2->Add(h_fE_mult2,-1);

rnpo->fChain->Draw("p_Etot>>h_pE_mult3","p_seg>0 && p_clustMult==3");
rnpo->fChain->Draw("f_Etot>>h_fE_mult3","f_seg>0 && f_clustMult==3","same");
h_pE_mult3->Sumw2();
TH1F *h_RnE_mult3 = (TH1F*)h_pE_mult3->Clone();
h_RnE_mult3->Add(h_fE_mult3,-1);

rnpo->fChain->Draw("p_Etot>>h_pE_mult4","p_seg>0 && p_clustMult==4");
rnpo->fChain->Draw("f_Etot>>h_fE_mult4","f_seg>0 && f_clustMult==4","same");
h_pE_mult4->Sumw2();
TH1F *h_RnE_mult4 = (TH1F*)h_pE_mult4->Clone();
h_RnE_mult4->Add(h_fE_mult4,-1);

rnpo->fChain->Draw("p_Etot>>h_pE_mult5","p_seg>0 && p_clustMult==5");
rnpo->fChain->Draw("f_Etot>>h_fE_mult5","f_seg>0 && f_clustMult==5","same");
h_pE_mult5->Sumw2();
TH1F *h_RnE_mult5 = (TH1F*)h_pE_mult5->Clone();
h_RnE_mult5->Add(h_fE_mult5,-1);

new TCanvas;
h_RnE_mult2->SetLineColor(kRed);
h_RnE_mult2->SetMarkerColor(kRed);
h_RnE_mult3->SetLineColor(8);
h_RnE_mult3->SetMarkerColor(8);
h_RnE_mult4->SetLineColor(kMagenta);
h_RnE_mult4->SetMarkerColor(kMagenta);
h_RnE_mult5->SetLineColor(9);
h_RnE_mult5->SetMarkerColor(9);

h_RnE_mult1->Draw();
h_RnE_mult2->Draw("same");
h_RnE_mult3->Draw("same");
h_RnE_mult4->Draw("same");
h_RnE_mult5->Draw("same");


//=====================================================================================

TH1F *h_pEn_mult1 = new TH1F("h_pEn_mult1","",150,0.5,1.15);
TH1F *h_fEn_mult1 = new TH1F("h_fEn_mult1","",150,0.5,1.15);

TH1F *h_pEn_mult2 = new TH1F("h_pEn_mult2","",150,0.5,1.15);
TH1F *h_fEn_mult2 = new TH1F("h_fEn_mult2","",150,0.5,1.15);

TH1F *h_pEn_mult3 = new TH1F("h_pEn_mult3","",150,0.5,1.15);
TH1F *h_fEn_mult3 = new TH1F("h_fEn_mult3","",150,0.5,1.15);

TH1F *h_pEn_mult4 = new TH1F("h_pEn_mult4","",150,0.5,1.15);
TH1F *h_fEn_mult4 = new TH1F("h_fEn_mult4","",150,0.5,1.15);

TH1F *h_pEn_mult5 = new TH1F("h_pEn_mult5","",150,0.5,1.15);
TH1F *h_fEn_mult5 = new TH1F("h_fEn_mult5","",150,0.5,1.15);

new TCanvas;

rnpo->fChain->Draw("p_E>>h_pEn_mult1","p_seg>0 && p_clustMult==1");
rnpo->fChain->Draw("f_E>>h_fEn_mult1","f_seg>0 && f_clustMult==1","same");
h_pEn_mult1->Sumw2();
TH1F *h_RnEn_mult1 = (TH1F*)h_pEn_mult1->Clone();
h_RnEn_mult1->Add(h_fEn_mult1,-1);

rnpo->fChain->Draw("p_E>>h_pEn_mult2","p_seg>0 && p_clustMult==2");
rnpo->fChain->Draw("f_E>>h_fEn_mult2","f_seg>0 && f_clustMult==2","same");
h_pEn_mult2->Sumw2();
TH1F *h_RnEn_mult2 = (TH1F*)h_pEn_mult2->Clone();
h_RnEn_mult2->Add(h_fEn_mult2,-1);

rnpo->fChain->Draw("p_E>>h_pEn_mult3","p_seg>0 && p_clustMult==3");
rnpo->fChain->Draw("f_E>>h_fEn_mult3","f_seg>0 && f_clustMult==3","same");
h_pEn_mult3->Sumw2();
TH1F *h_RnEn_mult3 = (TH1F*)h_pEn_mult3->Clone();
h_RnEn_mult3->Add(h_fEn_mult3,-1);

rnpo->fChain->Draw("p_E>>h_pEn_mult4","p_seg>0 && p_clustMult==4");
rnpo->fChain->Draw("f_E>>h_fEn_mult4","f_seg>0 && f_clustMult==4","same");
h_pEn_mult4->Sumw2();
TH1F *h_RnEn_mult4 = (TH1F*)h_pEn_mult4->Clone();
h_RnEn_mult4->Add(h_fEn_mult4,-1);

rnpo->fChain->Draw("p_E>>h_pEn_mult5","p_seg>0 && p_clustMult==5");
rnpo->fChain->Draw("f_E>>h_fEn_mult5","f_seg>0 && f_clustMult==5","same");
h_pEn_mult5->Sumw2();
TH1F *h_RnEn_mult5 = (TH1F*)h_pEn_mult5->Clone();
h_RnEn_mult5->Add(h_fEn_mult5,-1);

new TCanvas;
h_RnEn_mult2->SetLineColor(kRed);
h_RnEn_mult2->SetMarkerColor(kRed);
h_RnEn_mult3->SetLineColor(8);
h_RnEn_mult3->SetMarkerColor(8);
h_RnEn_mult4->SetLineColor(kMagenta);
h_RnEn_mult4->SetMarkerColor(kMagenta);
h_RnEn_mult5->SetLineColor(9);
h_RnEn_mult5->SetMarkerColor(9);

h_RnEn_mult1->Draw();
h_RnEn_mult2->Draw("same");
h_RnEn_mult3->Draw("same");
h_RnEn_mult4->Draw("same");
h_RnEn_mult5->Draw("same");



}

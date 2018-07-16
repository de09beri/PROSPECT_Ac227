//This macro will create histograms for Ac227 coincidences
//according to cell number

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
const int NUMCELLS = 154;
const double POLIFETIME = 2.57;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]

const int NUMEXCLUDECELLS = 31;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};

//const int NUMEXCLUDECELLS = 63;
//int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,21,23,24,27,28,32,34,40,41,42,44,52,55,56,68,69,70,79,83,84,86,97,98,102,111,112,115,122,125,126,127,130,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153};

void RnPoVsTime(double p_lowPSD, double d_lowPSD, double p_lowE, double d_lowE, double zLow, double zHigh, double timeBin, int dtFit){

	const double TIMEBREAK = timeBin*(3.6e6);	//[ms]

	//Set up bins and ranges for histograms
	double dtMin = 0.0, dtMax = TIMEWINDOW;	//[ms]
	int numDtBins = (dtMax - dtMin)/0.1;
	
	double PSDMin = 0.15, PSDMax = 0.37;
	int numPSDBins = 150; 

	double EnMin = 0.49, EnMax = 1.16;		//[MeVee]
	int numEnBins = (EnMax - EnMin)/0.005;	// 5 keV bins
	
	double dzMin = -250, dzMax = 250;	//[mm]
	int numDzBins = (dzMax - dzMin)/2.5;	//0.25 cm bins

	double posMin = -800, posMax = 800;	//[mm]
	int numPosBins = (posMax - posMin)/10;	//1 cm bins 

	TFile *histFile = new TFile(Form("%s/Ac227_HistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	TH1F *hSelectDt, 		*hBGDt,        *hRnPoDt;
	TH1F *hSelectPromptPSD, *hBGPromptPSD, *hRnPSD,	 *hSelectDelayPSD, *hBGDelayPSD, *hPoPSD;
	TH1F *hSelectPromptEn,  *hBGPromptEn,  *hRnEn,   *hSelectDelayEn,  *hBGDelayEn,  *hPoEn;
	TH1F *hSelectPromptPos, *hBGPromptPos, *hRnPos,  *hSelectDelayPos, *hBGDelayPos, *hPoPos;
	TH1F *hSelectDz,        *hBGDz,        *hRnPoDz;

	TH1F *hSelectPromptTotEn, *hBGPromptTotEn, *hRnTotEn;

	TH2F *hSelectPSDvsEn, 			*hBGPSDvsEn, 		   *hRnPoPSDvsEn;
	TH2F *hSelectDelayEnvsPromptEn, *hBGDelayEnvsPromptEn, *hPoEnvsRnEn;

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnCB,    *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	vector<double> vRate,   vRateErr;
	vector<double> vTotEff, vTotEffErr;
	vector<double> vLifetime,   vLifetimeErr;
	vector<double> vLivetime,   vTimestamp;
	vector<double> vPoPSDMean,  vPoPSDMeanErr,  vPoPSDSigma,  vPoPSDSigmaErr;
	vector<double> vPoEnMean,   vPoEnMeanErr,   vPoEnSigma,   vPoEnSigmaErr;
	vector<double> vPoPosMean,  vPoPosMeanErr,  vPoPosSigma,  vPoPosSigmaErr;
	vector<double> vRnPoDzMean, vRnPoDzMeanErr, vRnPoDzSigma, vRnPoDzSigmaErr;

	//---------------------------------------------------------------------------------
	RNPO *rnpo = new RNPO();
	
	bool exclude;
	Long64_t numEntries = Long64_t(rnpo->fChain->GetEntries());
	printf("Number of Ac-227 candidates: %lld \n",numEntries);

	double promptLowPSDCut, promptHighPSDCut, promptLowEnCut, promptHighEnCut;
	double delayLowPSDCut,  delayHighPSDCut,  delayLowEnCut,  delayHighEnCut;
	double dzCut;

	// Get cut values
	rnpo->GetEntry(0);
    promptLowPSDCut  = (p_lowPSD < rnpo->p_PSDCut[0]) ? p_lowPSD : rnpo->p_PSDCut[0];
    promptHighPSDCut = rnpo->p_PSDCut[1];
    delayLowPSDCut   = (d_lowPSD< rnpo->d_PSDCut[0]) ? d_lowPSD : rnpo->d_PSDCut[0];
    delayHighPSDCut  = rnpo->d_PSDCut[1];
    promptLowEnCut   = (p_lowE < rnpo->p_ECut[0]) ? p_lowE : rnpo->p_ECut[0];
    promptHighEnCut  = rnpo->p_ECut[1];
    delayLowEnCut    = (d_lowE < rnpo->d_ECut[0]) ? d_lowE : rnpo->d_ECut[0];
    delayHighEnCut   = rnpo->d_ECut[1];
    dzCut            = rnpo->dzCut;

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Loop over all events, splitting into time periods of 24 hours

	Long64_t IDX = 0;
	int numTimeBin = 0;
	while(IDX<numEntries){

		//---------------------------------------------------------------------------------
		//Initialize histograms
		hSelectDt 			= new TH1F(Form("hSelectDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
		hBGDt 				= new TH1F(Form("hBGDt_%i",numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
	
		hSelectPromptPSD 	= new TH1F(Form("hSelectPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGPromptPSD 		= new TH1F(Form("hBGPromptPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectDelayPSD 	= new TH1F(Form("hSelectDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
		hBGDelayPSD 		= new TH1F(Form("hBGDelayPSD_%i",numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

		hSelectPromptEn 	= new TH1F(Form("hSelectPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptEn 		= new TH1F(Form("hBGPromptEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
	
		hSelectDelayEn 		= new TH1F(Form("hSelectDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGDelayEn 			= new TH1F(Form("hBGDelayEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPromptPos 	= new TH1F(Form("hSelectPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGPromptPos 		= new TH1F(Form("hBGPromptPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDelayPos		= new TH1F(Form("hSelectDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
		hBGDelayPos 		= new TH1F(Form("hBGDelayPos_%i",numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

		hSelectDz 			= new TH1F(Form("hSelectDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
		hBGDz 				= new TH1F(Form("hBGDz_%i",numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

		hSelectPromptTotEn 	= new TH1F(Form("hSelectPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
		hBGPromptTotEn 		= new TH1F(Form("hBGPromptTotEn_%i",numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

		hSelectPSDvsEn 		= new TH2F(Form("hSelectPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
		hBGPSDvsEn 			= new TH2F(Form("hBGPSDvsEn_%i",numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

		hSelectDelayEnvsPromptEn 	= new TH2F(Form("hSelectDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		hBGDelayEnvsPromptEn 		= new TH2F(Form("hBGDelayEnvsPromptEn_%i",numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);

		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

		double lastOCSTime = 0.0, OCSTime = 0.0;
		double lastMuonVetoTime = 0.0, muonVetoTime = 0.0;

		for(Long64_t i=IDX;i<numEntries;i++){
			if(i%100000==0) printf("Event: %lld \n",i);
			rnpo->GetEntry(i);

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate livetime and weighted timestamp
			if(rnpo->d_t < lastTime){ 
				livetime += lastTime*(1e-6);		//livetime in ms	
				livetime = livetime - (lastOCSTime*(1e-6));
//				livetime = livetime - (lastMuonVetoTime*(1e-6));
	
				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;
			}

			if(livetime>TIMEBREAK){
				vLivetime.push_back(livetime);

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	

				IDX = i;
				break;
			}

			lastTime = rnpo->d_t;
			lastOCSTime = rnpo->OCSVeto_t;
//			lastMuonVetoTime = rnpo->muonVeto_t;
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);
				livetime = livetime - (lastOCSTime*(1e-6));
//				livetime = livetime - (lastMuonVetoTime*(1e-6));

				sumWeightedTimestamp += lastRunTime * ((lastRunTime/2.0)+lastTimestamp);
				sumRunTime += lastRunTime;

				vLivetime.push_back(livetime);

				tstamp = sumWeightedTimestamp/sumRunTime;
				vTimestamp.push_back(tstamp);	

				IDX = i+1;
			}

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Fill histograms

			seg = rnpo->d_seg;
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
			if(exclude) continue;

			if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;	
			if(rnpo->d_z < zLow || rnpo->d_z > zHigh) continue;

			if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut && rnpo->p_z>zLow && rnpo->p_z<zHigh){	
				dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
				dz = rnpo->d_z - rnpo->p_z;

				hSelectDt->Fill(dt);
				hSelectPromptPSD->Fill(rnpo->p_PSD);
				hSelectDelayPSD->Fill(rnpo->d_PSD);
				hSelectPromptEn->Fill(rnpo->p_E);
				hSelectDelayEn->Fill(rnpo->d_E);
				hSelectPromptTotEn->Fill(rnpo->p_Etot);	
				hSelectPromptPos->Fill(rnpo->p_z);
				hSelectDelayPos->Fill(rnpo->d_z);
				hSelectDz->Fill(dz);

				hSelectPSDvsEn->Fill(rnpo->p_E,rnpo->p_PSD);
				hSelectPSDvsEn->Fill(rnpo->d_E,rnpo->d_PSD);
				hSelectDelayEnvsPromptEn->Fill(rnpo->p_E,rnpo->d_E);		
			}
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut && rnpo->f_z>zLow && rnpo->f_z<zHigh){	
				dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
				dz = rnpo->d_z - rnpo->f_z;

				hBGDt->Fill(dt);
				hBGPromptPSD->Fill(rnpo->f_PSD);
				hBGDelayPSD->Fill(rnpo->d_PSD);
				hBGPromptEn->Fill(rnpo->f_E);
				hBGDelayEn->Fill(rnpo->d_E);
				hBGPromptTotEn->Fill(rnpo->f_Etot);	
				hBGPromptPos->Fill(rnpo->f_z);
				hBGDelayPos->Fill(rnpo->d_z);
				hBGDz->Fill(dz);

				hBGPSDvsEn->Fill(rnpo->f_E,rnpo->f_PSD);
				hBGPSDvsEn->Fill(rnpo->d_E,rnpo->d_PSD);
				hBGDelayEnvsPromptEn->Fill(rnpo->f_E,rnpo->d_E);		
			}

		}	//end for loop over events

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

		hSelectDt->Sumw2();
		hRnPoDt = (TH1F*)hSelectDt->Clone();
		hRnPoDt->SetName(Form("hRnPoDt_%i",numTimeBin));
		hBGDt->Sumw2();
		hRnPoDt->Add(hBGDt,-1);

		hRnPSD = (TH1F*)hSelectPromptPSD->Clone();
		hRnPSD->SetName(Form("hRnPSD_%i",numTimeBin));
		hRnPSD->Sumw2();
		hRnPSD->Add(hBGPromptPSD,-1);

		hPoPSD = (TH1F*)hSelectDelayPSD->Clone();
		hPoPSD->SetName(Form("hPoPSD_%i",numTimeBin));
		hPoPSD->Sumw2();
		hPoPSD->Add(hBGDelayPSD,-1);

		hRnEn = (TH1F*)hSelectPromptEn->Clone();
		hRnEn->SetName(Form("hRnEn_%i",numTimeBin));
		hRnEn->Sumw2();
		hRnEn->Add(hBGPromptEn,-1);

		hRnTotEn = (TH1F*)hSelectPromptTotEn->Clone();
		hRnTotEn->SetName(Form("hRnTotEn_%i",numTimeBin));
		hRnTotEn->Sumw2();
		hRnTotEn->Add(hBGPromptTotEn,-1);

		hPoEn = (TH1F*)hSelectDelayEn->Clone();
		hPoEn->SetName(Form("hPoEn_%i",numTimeBin));
		hPoEn->Sumw2();
		hPoEn->Add(hBGDelayEn,-1);
	
		hRnPos = (TH1F*)hSelectPromptPos->Clone();
		hRnPos->SetName(Form("hRnPos_%i",numTimeBin));
		hRnPos->Sumw2();
		hRnPos->Add(hBGPromptPos,-1);

		hPoPos = (TH1F*)hSelectDelayPos->Clone();
		hPoPos->SetName(Form("hPoPos_%i",numTimeBin));
		hPoPos->Sumw2();
		hPoPos->Add(hBGDelayPos,-1);

		hRnPoDz = (TH1F*)hSelectDz->Clone();
		hRnPoDz->SetName(Form("hRnPoDz_%i",numTimeBin));
		hRnPoDz->Sumw2();
		hRnPoDz->Add(hBGDz,-1);

		hRnPoPSDvsEn = (TH2F*)hSelectPSDvsEn->Clone();
		hRnPoPSDvsEn->SetName(Form("hRnPoPSDvsEn_%i",numTimeBin));
		hRnPoPSDvsEn->Add(hBGPSDvsEn,-1);	

		hPoEnvsRnEn = (TH2F*)hSelectDelayEnvsPromptEn->Clone();
		hPoEnvsRnEn->SetName(Form("hPoEnvsRnEn_%i",numTimeBin));
		hPoEnvsRnEn->Add(hBGDelayEnvsPromptEn,-1);

		//---------------------------------------------------------------------------------
		//Calculate results
		printf("=============== Calculating Results =============== \n"); 

		double dtBinWidth = (dtMax - dtMin)/(double)numDtBins;

		double promptPSDEff, delayPSDEff, promptPSDEffErr, delayPSDEffErr;
		double promptEnEff,  delayEnEff,  promptEnEffErr,  delayEnEffErr;
		double dzEff, dzEffErr;
		double totEff, totEffErr;

		double NAlpha, NAlphaErr, lifetime, lifetimeErr;
		double rate, rateErr;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Fit distributions
		fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),dtBinWidth*dtFit,dtMax);
		fRnPoDtExp->SetParameter(1,POLIFETIME);
		hRnPoDt->Fit(fRnPoDtExp,"RQ0");

		double fitPSDMin = 0.20;
		fRnPSDGaus = new TF1("fRnPSDGaus","gaus",fitPSDMin,PSDMax);
		hRnPSD->Fit(fRnPSDGaus,"RQ0");
		fRnPSDGaus->SetRange(PSDMin,PSDMax);

		fPoPSDGaus = new TF1("fPoPSDGaus","gaus",fitPSDMin,PSDMax);
		hPoPSD->Fit(fPoPSDGaus,"RQ0");
		fPoPSDGaus->SetRange(PSDMin,PSDMax);
	
		double fitRnEnMin = 0.57;	
		fRnEnCB = new TF1("fRnEnCB","crystalball",fitRnEnMin,EnMax);
		fRnEnCB->SetParameter(1,0.7);
		fRnEnCB->SetParameter(2,0.04);
		fRnEnCB->SetParameter(3,-1.7);
		fRnEnCB->SetParameter(4,1.2);
		hRnEn->Fit(fRnEnCB,"RQ0");
		fRnEnCB->SetRange(EnMin,EnMax);

		double fitPoEnMin = 0.7;
		fPoEnGaus = new TF1("fPoEnGaus","gaus",fitPoEnMin,EnMax);
		hPoEn->Fit(fPoEnGaus,"RQ0");
		fPoEnGaus->SetRange(EnMin,EnMax);

		fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",dzMin,dzMax);
		hRnPoDz->Fit(fRnPoDzGaus,"RQ0");

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate efficiencies
		promptPSDEff = fRnPSDGaus->Integral(promptLowPSDCut,promptHighPSDCut)/fRnPSDGaus->Integral(PSDMin,PSDMax);
		promptPSDEffErr = sqrt((promptPSDEff*(1-promptPSDEff))/hRnPSD->GetEntries()); 
			
		delayPSDEff = fPoPSDGaus->Integral(delayLowPSDCut,delayHighPSDCut)/fPoPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1-delayPSDEff))/hPoPSD->GetEntries()); 
		
		promptEnEff = fRnEnCB->Integral(promptLowEnCut,promptHighEnCut)/fRnEnCB->Integral(EnMin,EnMax);
		promptEnEffErr = sqrt((promptEnEff*(1-promptEnEff))/hRnEn->GetEntries());

		delayEnEff = fPoEnGaus->Integral(delayLowEnCut,delayHighEnCut)/fPoEnGaus->Integral(EnMin,EnMax);
		delayEnEffErr = sqrt((delayEnEff*(1-delayEnEff))/hPoEn->GetEntries()); 

		dzEff = fRnPoDzGaus->Integral(-dzCut,dzCut)/fRnPoDzGaus->Integral(dzMin,dzMax);
		dzEffErr = sqrt((dzEff*(1-dzEff))/hRnPoDz->GetEntries());

		
		totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
		totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//Calculate rate
		NAlpha = fRnPoDtExp->GetParameter(0);
		NAlphaErr = fRnPoDtExp->GetParError(0);
		lifetime = fRnPoDtExp->GetParameter(1);
		lifetimeErr = fRnPoDtExp->GetParError(1);

		rate = (NAlpha/(livetime*totEff))*(1e3);		//Hz
		rateErr = rate * sqrt(pow(NAlphaErr/NAlpha,2) + pow(totEffErr/totEff,2));

		printf("Rate: %.4f +/- %.4f \n",rate,rateErr);

		//---------------------------------------------------------------------------------
		//Populate vectors
		vRate.push_back(rate);
		vRateErr.push_back(rateErr);
		vTotEff.push_back(totEff);
		vTotEffErr.push_back(totEffErr);
		vLifetime.push_back(lifetime);
		vLifetimeErr.push_back(lifetimeErr);
		
		vPoPSDMean.push_back(fPoPSDGaus->GetParameter(1));	
		vPoPSDMeanErr.push_back(fPoPSDGaus->GetParError(1));
		vPoPSDSigma.push_back(fPoPSDGaus->GetParameter(2));
		vPoPSDSigmaErr.push_back(fPoPSDGaus->GetParError(2));

		vPoEnMean.push_back(fPoEnGaus->GetParameter(1));
		vPoEnMeanErr.push_back(fPoEnGaus->GetParError(1));
		vPoEnSigma.push_back(fPoEnGaus->GetParameter(2));
		vPoEnSigmaErr.push_back(fPoEnGaus->GetParError(2));

		vPoPosMean.push_back(hPoPos->GetMean());
		vPoPosMeanErr.push_back(hPoPos->GetMeanError());
		vPoPosSigma.push_back(hPoPos->GetRMS());
		vPoPosSigmaErr.push_back(hPoPos->GetRMSError());

		vRnPoDzMean.push_back(fRnPoDzGaus->GetParameter(1));
		vRnPoDzMeanErr.push_back(fRnPoDzGaus->GetParError(1));
		vRnPoDzSigma.push_back(fRnPoDzGaus->GetParameter(2));
		vRnPoDzSigmaErr.push_back(fRnPoDzGaus->GetParError(2));

		numTimeBin++;
	}	//end while loop IDX < numEntries

	histFile->Write();
	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numPt = vRate.size();
	double x[numPt], xErr[numPt];
	double y[1], yErr[1];

	TGraphErrors *grRate 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grTotEff 		= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grLifetime 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPSDSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean   	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean 	= new TGraphErrors(numPt,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma = new TGraphErrors(numPt,x,y,xErr,yErr);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double time = vTimestamp[i];

		grRate->SetPoint(i,time,vRate[i]);
		grRate->SetPointError(i,0,vRateErr[i]);

		grTotEff->SetPoint(i,time,vTotEff[i]);
		grTotEff->SetPointError(i,0,vTotEffErr[i]);

		grLifetime->SetPoint(i,time,vLifetime[i]);
		grLifetime->SetPointError(i,0,vLifetimeErr[i]);

		grPoPSDMean->SetPoint(i,time,vPoPSDMean[i]);
		grPoPSDMean->SetPointError(i,0,vPoPSDMeanErr[i]);

		grPoPSDSigma->SetPoint(i,time,vPoPSDSigma[i]);
		grPoPSDSigma->SetPointError(i,0,vPoPSDSigmaErr[i]);

		grPoEnMean->SetPoint(i,time,vPoEnMean[i]);
		grPoEnMean->SetPointError(i,0,vPoEnMeanErr[i]);
	
		grPoEnSigma->SetPoint(i,time,vPoEnSigma[i]);
		grPoEnSigma->SetPointError(i,0,vPoEnSigmaErr[i]);

		grPoPosMean->SetPoint(i,time,vPoPosMean[i]);
		grPoPosMean->SetPointError(i,0,vPoPosMeanErr[i]);

		grPoPosSigma->SetPoint(i,time,vPoPosSigma[i]);
		grPoPosSigma->SetPointError(i,0,vPoPosSigmaErr[i]);		

		grRnPoDzMean->SetPoint(i,time,vRnPoDzMean[i]);
		grRnPoDzMean->SetPointError(i,0,vRnPoDzMeanErr[i]);

		grRnPoDzSigma->SetPoint(i,time,vRnPoDzSigma[i]);
		grRnPoDzSigma->SetPointError(i,0,vRnPoDzSigmaErr[i]);
	}	//end for loop to populate TGraphs

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_GraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");
	
	grRate->Write("grRate");
	grTotEff->Write("grTotEff");	
	grLifetime->Write("grLifetime");
	grPoPSDMean->Write("grPoPSDMean");
	grPoPSDSigma->Write("grPoPSDSigma");
	grPoEnMean->Write("grPoEnMean");
	grPoEnSigma->Write("grPoEnSigma");
	grPoPosMean->Write("grPoPosMean");
	grPoPosSigma->Write("grPoPosSigma");
	grRnPoDzMean->Write("grRnPoDzMean");
	grRnPoDzSigma->Write("grRnPoDzSigma");	

	graphFile->Close();

}	//end void RnPoVsTime


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
const int NUMCOLS = 14;
const double POLIFETIME = 2.57;   //[ms]
const double TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;    //[ms]

const int NUMEXCLUDECELLS = 33;
int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,8,9,10,11,12,13,18,21,23,24,26,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};

void RnPoColVsTime(double p_lowPSD, double d_lowPSD, double p_lowE, double d_lowE, double timeBin, int dtFit){
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

	TFile *histFile = new TFile(Form("%s/Ac227_ColHistsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	//---------------------------------------------------------------------------------
	TH1F *hSelectDt[NUMCOLS], 		 *hBGDt[NUMCOLS],        *hRnPoDt[NUMCOLS];
	TH1F *hSelectPromptPSD[NUMCOLS], *hBGPromptPSD[NUMCOLS], *hRnPSD[NUMCOLS],	*hSelectDelayPSD[NUMCOLS], *hBGDelayPSD[NUMCOLS], *hPoPSD[NUMCOLS];
	TH1F *hSelectPromptEn[NUMCOLS],  *hBGPromptEn[NUMCOLS],  *hRnEn[NUMCOLS],   *hSelectDelayEn[NUMCOLS],  *hBGDelayEn[NUMCOLS],  *hPoEn[NUMCOLS];
	TH1F *hSelectPromptPos[NUMCOLS], *hBGPromptPos[NUMCOLS], *hRnPos[NUMCOLS],  *hSelectDelayPos[NUMCOLS], *hBGDelayPos[NUMCOLS], *hPoPos[NUMCOLS];
	TH1F *hSelectDz[NUMCOLS],        *hBGDz[NUMCOLS],        *hRnPoDz[NUMCOLS];

	TH1F *hSelectPromptTotEn[NUMCOLS], *hBGPromptTotEn[NUMCOLS], *hRnTotEn[NUMCOLS];

	TH2F *hSelectPSDvsEn[NUMCOLS], 			 *hBGPSDvsEn[NUMCOLS], 		     *hRnPoPSDvsEn[NUMCOLS];
	TH2F *hSelectDelayEnvsPromptEn[NUMCOLS], *hBGDelayEnvsPromptEn[NUMCOLS], *hPoEnvsRnEn[NUMCOLS];

	//---------------------------------------------------------------------------------
	TF1 *fRnPoDtExp;
	TF1 *fRnPSDGaus, *fPoPSDGaus;
	TF1 *fRnEnCB,    *fPoEnGaus;
	TF1 *fRnPoDzGaus;

	//---------------------------------------------------------------------------------
	int ncol = NUMCOLS, nt = 500;

	double aRate[ncol][nt], aRateErr[ncol][nt];
	double aTotEff[ncol][nt], aTotEffErr[ncol][nt];
	double aLifetime[ncol][nt],   aLifetimeErr[ncol][nt];
	double aPoPSDMean[ncol][nt],  aPoPSDMeanErr[ncol][nt],  aPoPSDSigma[ncol][nt],  aPoPSDSigmaErr[ncol][nt];
	double aPoEnMean[ncol][nt],   aPoEnMeanErr[ncol][nt],   aPoEnSigma[ncol][nt],   aPoEnSigmaErr[ncol][nt];
	double aPoPosMean[ncol][nt],  aPoPosMeanErr[ncol][nt],  aPoPosSigma[ncol][nt],  aPoPosSigmaErr[ncol][nt];
	double aRnPoDzMean[ncol][nt], aRnPoDzMeanErr[ncol][nt], aRnPoDzSigma[ncol][nt], aRnPoDzSigmaErr[ncol][nt];

	vector<double> vLivetime,   vTimestamp;

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
		for(int i=0;i<NUMCOLS;i++){
			hSelectDt[i] 			= new TH1F(Form("hSelectDt_Col%i_%i",i,numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
			hBGDt[i] 				= new TH1F(Form("hBGDt_Col%i_%i",i,numTimeBin),";dt [ms];Counts/0.1 ms",numDtBins,dtMin,dtMax);
		
			hSelectPromptPSD[i] 	= new TH1F(Form("hSelectPromptPSD_Col%i_%i",i,numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
			hBGPromptPSD[i] 		= new TH1F(Form("hBGPromptPSD_Col%i_%i",i,numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

			hSelectDelayPSD[i] 		= new TH1F(Form("hSelectDelayPSD_Col%i_%i",i,numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	
			hBGDelayPSD[i] 			= new TH1F(Form("hBGDelayPSD_Col%i_%i",i,numTimeBin),";PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);	

			hSelectPromptEn[i] 		= new TH1F(Form("hSelectPromptEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
			hBGPromptEn[i] 			= new TH1F(Form("hBGPromptEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
			
			hSelectDelayEn[i] 		= new TH1F(Form("hSelectDelayEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
			hBGDelayEn[i] 			= new TH1F(Form("hBGDelayEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

			hSelectPromptPos[i] 	= new TH1F(Form("hSelectPromptPos_Col%i_%i",i,numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
			hBGPromptPos[i] 		= new TH1F(Form("hBGPromptPos_Col%i_%i",i,numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

			hSelectDelayPos[i]		= new TH1F(Form("hSelectDelayPos_Col%i_%i",i,numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);
			hBGDelayPos[i] 			= new TH1F(Form("hBGDelayPos_Col%i_%i",i,numTimeBin),";z [mm];Counts/cm",numPosBins,posMin,posMax);

			hSelectDz[i] 			= new TH1F(Form("hSelectDz_Col%i_%i",i,numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);
			hBGDz[i] 				= new TH1F(Form("hBGDz_Col%i_%i",i,numTimeBin),";z_{Po} - z_{Rn} [mm];Counts/0.25 cm",numDzBins,dzMin,dzMax);

			hSelectPromptTotEn[i] 	= new TH1F(Form("hSelectPromptTotEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	
			hBGPromptTotEn[i] 		= new TH1F(Form("hBGPromptTotEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];Counts/5 keV",numEnBins,EnMin,EnMax);	

			hSelectPSDvsEn[i] 		= new TH2F(Form("hSelectPSDvsEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);
			hBGPSDvsEn[i] 			= new TH2F(Form("hBGPSDvsEn_Col%i_%i",i,numTimeBin),";Energy [MeVee];PSD [arb]",numEnBins,EnMin,EnMax,numPSDBins,PSDMin,PSDMax);

			hSelectDelayEnvsPromptEn[i] 	= new TH2F(Form("hSelectDelayEnvsPromptEn_Col%i_%i",i,numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
			hBGDelayEnvsPromptEn[i] 		= new TH2F(Form("hBGDelayEnvsPromptEn_Col%i_%i",i,numTimeBin),";Rn Energy [MeVee];Po Energy [MeVee]",numEnBins,EnMin,EnMax,numEnBins,EnMin,EnMax);
		}

		//---------------------------------------------------------------------------------
		//Fill histograms
		printf("=============== Filling Histograms =============== \n"); 
	
		int seg, col;
		double dt, dz;
		
		double livetime = 0.0, tstamp;
		double lastTime = 0.0, lastRunTime = 0.0, lastTimestamp = 0.0;	
		double sumWeightedTimestamp = 0.0, sumRunTime = 0.0;

		for(Long64_t i=IDX;i<numEntries;i++){
			if(i%100000==0) printf("Event: %lld \n",i);
			rnpo->GetEntry(i);

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate livetime and weighted timestamp
			if(rnpo->d_t < lastTime){ 
				livetime += lastTime*(1e-6);		//livetime in ms	
	
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
			lastRunTime = ((TVectorD*)rnpo->fChain->GetCurrentFile()->Get("runtime"))->Norm1();		//seconds
			lastTimestamp = rnpo->tstamp;			//epoch seconds	

			//if we are at the last entry	
			if(i == (numEntries-1)){ 
				livetime += lastTime*(1e-6);

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
			col = seg%NUMCOLS;
			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), seg) != end(ExcludeCellArr);
			if(exclude) continue;

			if(rnpo->d_PSD < delayLowPSDCut || rnpo->d_E < delayLowEnCut) continue;	

			if(rnpo->p_seg > -1 && rnpo->p_PSD>promptLowPSDCut && rnpo->p_E>promptLowEnCut){	
				dt = (rnpo->d_t - rnpo->p_t)*(1e-6);	//convert ns to ms	
				dz = rnpo->d_z - rnpo->p_z;

				hSelectDt[col]->Fill(dt);
				hSelectPromptPSD[col]->Fill(rnpo->p_PSD);
				hSelectDelayPSD[col]->Fill(rnpo->d_PSD);
				hSelectPromptEn[col]->Fill(rnpo->p_E);
				hSelectDelayEn[col]->Fill(rnpo->d_E);
				hSelectPromptTotEn[col]->Fill(rnpo->p_Etot);	
				hSelectPromptPos[col]->Fill(rnpo->p_z);
				hSelectDelayPos[col]->Fill(rnpo->d_z);
				hSelectDz[col]->Fill(dz);

				hSelectPSDvsEn[col]->Fill(rnpo->p_E,rnpo->p_PSD);
				hSelectPSDvsEn[col]->Fill(rnpo->d_E,rnpo->d_PSD);
				hSelectDelayEnvsPromptEn[col]->Fill(rnpo->p_E,rnpo->d_E);		
			}
			if(rnpo->f_seg > -1 && rnpo->f_PSD>promptLowPSDCut && rnpo->f_E>promptLowEnCut){	
				dt = (rnpo->f_t - rnpo->d_t)*(1e-6) - TIMEOFFSET;	
				dz = rnpo->d_z - rnpo->f_z;

				hBGDt[col]->Fill(dt);
				hBGPromptPSD[col]->Fill(rnpo->f_PSD);
				hBGDelayPSD[col]->Fill(rnpo->d_PSD);
				hBGPromptEn[col]->Fill(rnpo->f_E);
				hBGDelayEn[col]->Fill(rnpo->d_E);
				hBGPromptTotEn[col]->Fill(rnpo->f_Etot);	
				hBGPromptPos[col]->Fill(rnpo->f_z);
				hBGDelayPos[col]->Fill(rnpo->d_z);
				hBGDz[col]->Fill(dz);

				hBGPSDvsEn[col]->Fill(rnpo->f_E,rnpo->f_PSD);
				hBGPSDvsEn[col]->Fill(rnpo->d_E,rnpo->d_PSD);
				hBGDelayEnvsPromptEn[col]->Fill(rnpo->f_E,rnpo->d_E);		
			}

		}	//end for loop over events

		printf("Time bin: %i  |  Livetime: %f hrs \n",numTimeBin,livetime*(2.778e-7));

		//---------------------------------------------------------------------------------
		//Subtract histograms
		printf("=============== Subtracting Histograms =============== \n"); 

		for(int i=0;i<NUMCOLS;i++){
			hRnPoDt[i] = (TH1F*)hSelectDt[i]->Clone();
			hRnPoDt[i]->SetName(Form("hRnPoDt_Col%i_%i",i,numTimeBin));
			hRnPoDt[i]->Sumw2();
			hRnPoDt[i]->Add(hBGDt[i],-1);

			hRnPSD[i] = (TH1F*)hSelectPromptPSD[i]->Clone();
			hRnPSD[i]->SetName(Form("hRnPSD_Col%i_%i",i,numTimeBin));
			hRnPSD[i]->Sumw2();
			hRnPSD[i]->Add(hBGPromptPSD[i],-1);

			hPoPSD[i] = (TH1F*)hSelectDelayPSD[i]->Clone();
			hPoPSD[i]->SetName(Form("hPoPSD_Col%i_%i",i,numTimeBin));
			hPoPSD[i]->Sumw2();
			hPoPSD[i]->Add(hBGDelayPSD[i],-1);

			hRnEn[i] = (TH1F*)hSelectPromptEn[i]->Clone();
			hRnEn[i]->SetName(Form("hRnEn_Col%i_%i",i,numTimeBin));
			hRnEn[i]->Sumw2();
			hRnEn[i]->Add(hBGPromptEn[i],-1);

			hRnTotEn[i] = (TH1F*)hSelectPromptTotEn[i]->Clone();
			hRnTotEn[i]->SetName(Form("hRnTotEn_Col%i_%i",i,numTimeBin));
			hRnTotEn[i]->Sumw2();
			hRnTotEn[i]->Add(hBGPromptTotEn[i],-1);

			hPoEn[i] = (TH1F*)hSelectDelayEn[i]->Clone();
			hPoEn[i]->SetName(Form("hPoEn_Col%i_%i",i,numTimeBin));
			hPoEn[i]->Sumw2();
			hPoEn[i]->Add(hBGDelayEn[i],-1);
		
			hRnPos[i] = (TH1F*)hSelectPromptPos[i]->Clone();
			hRnPos[i]->SetName(Form("hRnPos_Col%i_%i",i,numTimeBin));
			hRnPos[i]->Sumw2();
			hRnPos[i]->Add(hBGPromptPos[i],-1);

			hPoPos[i] = (TH1F*)hSelectDelayPos[i]->Clone();
			hPoPos[i]->SetName(Form("hPoPos_Col%i_%i",i,numTimeBin));
			hPoPos[i]->Sumw2();
			hPoPos[i]->Add(hBGDelayPos[i],-1);

			hRnPoDz[i] = (TH1F*)hSelectDz[i]->Clone();
			hRnPoDz[i]->SetName(Form("hRnPoDz_Col%i_%i",i,numTimeBin));
			hRnPoDz[i]->Sumw2();
			hRnPoDz[i]->Add(hBGDz[i],-1);

			hRnPoPSDvsEn[i] = (TH2F*)hSelectPSDvsEn[i]->Clone();
			hRnPoPSDvsEn[i]->SetName(Form("hRnPoPSDvsEn_Col%i_%i",i,numTimeBin));
			hRnPoPSDvsEn[i]->Add(hBGPSDvsEn[i],-1);	

			hPoEnvsRnEn[i] = (TH2F*)hSelectDelayEnvsPromptEn[i]->Clone();
			hPoEnvsRnEn[i]->SetName(Form("hPoEnvsRnEn_Col%i_%i",i,numTimeBin));
			hPoEnvsRnEn[i]->Add(hBGDelayEnvsPromptEn[i],-1);
		}
		
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
		for(int i=0;i<NUMCOLS;i++){
			fRnPoDtExp = new TF1("fRnPoDtExp",Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth),dtBinWidth*dtFit,dtMax);
			fRnPoDtExp->SetParameter(1,POLIFETIME);
			hRnPoDt[i]->Fit(fRnPoDtExp,"RQ0");

			fRnPSDGaus = new TF1("fRnPSDGaus","gaus",PSDMin,PSDMax);
			hRnPSD[i]->Fit(fRnPSDGaus,"RQ0");

			fPoPSDGaus = new TF1("fPoPSDGaus","gaus",PSDMin,PSDMax);
			hPoPSD[i]->Fit(fPoPSDGaus,"RQ0");
			
			fRnEnCB = new TF1("fRnEnCB","crystalball",EnMin,EnMax);
			fRnEnCB->SetParameter(1,0.7);
			fRnEnCB->SetParameter(2,0.04);
			fRnEnCB->SetParameter(3,-1.7);
			fRnEnCB->SetParameter(4,1.2);
			hRnEn[i]->Fit(fRnEnCB,"RQ0");

			fPoEnGaus = new TF1("fPoEnGaus","gaus",EnMin,EnMax);
			hPoEn[i]->Fit(fPoEnGaus,"RQ0");

			fRnPoDzGaus = new TF1("fRnPoDzGaus","gaus",dzMin,dzMax);
			hRnPoDz[i]->Fit(fRnPoDzGaus,"RQ0");

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate efficiencies
			promptPSDEff = fRnPSDGaus->Integral(promptLowPSDCut,promptHighPSDCut)/fRnPSDGaus->Integral(PSDMin,PSDMax);
			promptPSDEffErr = sqrt((promptPSDEff*(1-promptPSDEff))/hRnPSD[i]->GetEntries()); 
				
			delayPSDEff = fPoPSDGaus->Integral(delayLowPSDCut,delayHighPSDCut)/fPoPSDGaus->Integral(PSDMin,PSDMax);
			delayPSDEffErr = sqrt((delayPSDEff*(1-delayPSDEff))/hPoPSD[i]->GetEntries()); 
			
			promptEnEff = fRnEnCB->Integral(promptLowEnCut,promptHighEnCut)/fRnEnCB->Integral(EnMin,EnMax);
			promptEnEffErr = sqrt((promptEnEff*(1-promptEnEff))/hRnEn[i]->GetEntries());

			delayEnEff = fPoEnGaus->Integral(delayLowEnCut,delayHighEnCut)/fPoEnGaus->Integral(EnMin,EnMax);
			delayEnEffErr = sqrt((delayEnEff*(1-delayEnEff))/hPoEn[i]->GetEntries()); 

			dzEff = fRnPoDzGaus->Integral(-dzCut,dzCut)/fRnPoDzGaus->Integral(dzMin,dzMax);
			dzEffErr = sqrt((dzEff*(1-dzEff))/hRnPoDz[i]->GetEntries());


			totEff = promptPSDEff * delayPSDEff * promptEnEff * delayEnEff * dzEff;	
			totEffErr = totEff * sqrt( pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(dzEffErr/dzEff,2) );

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			//Calculate rate
			NAlpha = fRnPoDtExp->GetParameter(0);
			NAlphaErr = fRnPoDtExp->GetParError(0);
			lifetime = fRnPoDtExp->GetParameter(1);
			lifetimeErr = fRnPoDtExp->GetParError(1);

			rate = (NAlpha/(livetime*totEff))*(1e6);		//mHz
			rateErr = rate * sqrt(pow(NAlphaErr/NAlpha,2) + pow(totEffErr/totEff,2));

			printf("Column %i Rate: %.4f +/- %.4f \n",i,rate,rateErr);

			//---------------------------------------------------------------------------------
			//Populate vectors
			aRate[i][numTimeBin] = rate;
			aRateErr[i][numTimeBin] = rateErr;
			aTotEff[i][numTimeBin] = totEff;
			aTotEffErr[i][numTimeBin] = totEffErr;
			aLifetime[i][numTimeBin] = lifetime;
			aLifetimeErr[i][numTimeBin] = lifetimeErr;
			
			aPoPSDMean[i][numTimeBin] = fPoPSDGaus->GetParameter(1);	
			aPoPSDMeanErr[i][numTimeBin] = fPoPSDGaus->GetParError(1);
			aPoPSDSigma[i][numTimeBin] = fPoPSDGaus->GetParameter(2);
			aPoPSDSigmaErr[i][numTimeBin] = fPoPSDGaus->GetParError(2);

			aPoEnMean[i][numTimeBin] = fPoEnGaus->GetParameter(1);
			aPoEnMeanErr[i][numTimeBin] = fPoEnGaus->GetParError(1);
			aPoEnSigma[i][numTimeBin] = fPoEnGaus->GetParameter(2);
			aPoEnSigmaErr[i][numTimeBin] = fPoEnGaus->GetParError(2);

			aPoPosMean[i][numTimeBin] = hPoPos[i]->GetMean();
			aPoPosMeanErr[i][numTimeBin] = hPoPos[i]->GetMeanError();
			aPoPosSigma[i][numTimeBin] = hPoPos[i]->GetRMS();
			aPoPosSigmaErr[i][numTimeBin] = hPoPos[i]->GetRMSError();

			aRnPoDzMean[i][numTimeBin] = fRnPoDzGaus->GetParameter(1);
			aRnPoDzMeanErr[i][numTimeBin] = fRnPoDzGaus->GetParError(1);
			aRnPoDzSigma[i][numTimeBin] = fRnPoDzGaus->GetParameter(2);
			aRnPoDzSigmaErr[i][numTimeBin] = fRnPoDzGaus->GetParError(2);
		}


		numTimeBin++;
	}	//end while loop IDX < numEntries

	histFile->Write();
	histFile->Close();

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	//---------------------------------------------------------------------------------
	//Initialize TGraphs for results

	int numPt = numTimeBin;
	double x[numPt], xErr[numPt];
	double y[1], yErr[1];

	TGraphErrors *grRate[NUMCOLS];
	TGraphErrors *grTotEff[NUMCOLS];
	TGraphErrors *grLifetime[NUMCOLS];
	TGraphErrors *grPoPSDMean[NUMCOLS],  *grPoPSDSigma[NUMCOLS];
	TGraphErrors *grPoEnMean[NUMCOLS],   *grPoEnSigma[NUMCOLS];
	TGraphErrors *grPoPosMean[NUMCOLS],  *grPoPosSigma[NUMCOLS];
	TGraphErrors *grRnPoDzMean[NUMCOLS], *grRnPoDzSigma[NUMCOLS];

	for(int i=0;i<NUMCOLS;i++){
		grRate[i] 		 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grTotEff[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grLifetime[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPSDMean[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPSDSigma[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoEnMean[i]    = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoEnSigma[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPosMean[i] 	 = new TGraphErrors(numPt,x,y,xErr,yErr);
		grPoPosSigma[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grRnPoDzMean[i]  = new TGraphErrors(numPt,x,y,xErr,yErr);
		grRnPoDzSigma[i] = new TGraphErrors(numPt,x,y,xErr,yErr);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//Fill TGraphs
	for(int i=0;i<numPt;i++){
		double time = vTimestamp[i];

		for(int j=0;j<NUMCOLS;j++){
			grRate[j]->SetPoint(i,time,aRate[j][i]);
			grRate[j]->SetPointError(i,0,aRateErr[j][i]);

			grTotEff[j]->SetPoint(i,time,aTotEff[j][i]);
			grTotEff[j]->SetPointError(i,0,aTotEffErr[j][i]);

			grLifetime[j]->SetPoint(i,time,aLifetime[j][i]);
			grLifetime[j]->SetPointError(i,0,aLifetimeErr[j][i]);

			grPoPSDMean[j]->SetPoint(i,time,aPoPSDMean[j][i]);
			grPoPSDMean[j]->SetPointError(i,0,aPoPSDMeanErr[j][i]);

			grPoPSDSigma[j]->SetPoint(i,time,aPoPSDSigma[j][i]);
			grPoPSDSigma[j]->SetPointError(i,0,aPoPSDSigmaErr[j][i]);

			grPoEnMean[j]->SetPoint(i,time,aPoEnMean[j][i]);
			grPoEnMean[j]->SetPointError(i,0,aPoEnMeanErr[j][i]);
		
			grPoEnSigma[j]->SetPoint(i,time,aPoEnSigma[j][i]);
			grPoEnSigma[j]->SetPointError(i,0,aPoEnSigmaErr[j][i]);

			grPoPosMean[j]->SetPoint(i,time,aPoPosMean[j][i]);
			grPoPosMean[j]->SetPointError(i,0,aPoPosMeanErr[j][i]);

			grPoPosSigma[j]->SetPoint(i,time,aPoPosSigma[j][i]);
			grPoPosSigma[j]->SetPointError(i,0,aPoPosSigmaErr[j][i]);		

			grRnPoDzMean[j]->SetPoint(i,time,aRnPoDzMean[j][i]);
			grRnPoDzMean[j]->SetPointError(i,0,aRnPoDzMeanErr[j][i]);

			grRnPoDzSigma[j]->SetPoint(i,time,aRnPoDzSigma[j][i]);
			grRnPoDzSigma[j]->SetPointError(i,0,aRnPoDzSigmaErr[j][i]);
		}
	}	//end for loop to populate TGraphs

	//---------------------------------------------------------------------------------
	//Write TGraphs to file
	TFile *graphFile = new TFile(Form("%s/Ac227_ColGraphsPerTime.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"RECREATE");

	for(int i=0;i<NUMCOLS;i++){
		grRate[i]->Write(Form("grRate_%i",i));
		grTotEff[i]->Write(Form("grTotEff_%i",i));	
		grLifetime[i]->Write(Form("grLifetime_%i",i));
		grPoPSDMean[i]->Write(Form("grPoPSDMean_%i",i));
		grPoPSDSigma[i]->Write(Form("grPoPSDSigma_%i",i));
		grPoEnMean[i]->Write(Form("grPoEnMean_%i",i));
		grPoEnSigma[i]->Write(Form("grPoEnSigma_%i",i));
		grPoPosMean[i]->Write(Form("grPoPosMean_%i",i));
		grPoPosSigma[i]->Write(Form("grPoPosSigma_%i",i));
		grRnPoDzMean[i]->Write(Form("grRnPoDzMean_%i",i));
		grRnPoDzSigma[i]->Write(Form("grRnPoDzSigma_%i",i));	
	}

	graphFile->Close();

}	//end void RnPoVsTime


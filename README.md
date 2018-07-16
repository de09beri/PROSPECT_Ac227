# Ac227Analysis
The Ac227 RnPo analysis for the PROSPECT AD happens in two steps.
The first step is the Ac227TreePlugin in the P2x framework. This plugin will create a TTree, TAc, containing information 
on candidate RnPo events.
The second step is analyzing these events to calculate rates and other interesting information.
This is done outside of the P2x framework with the macros described below.

## Creating RnPo Trees
The TAc trees are created using the Ac227TreePlugin. The energy, PSD, and dz cuts are defined in Ac227TreePlugin.hh. These are set fairly wide and should not have to be changed at this point. They can be tightened at a later step in the analysis.
	* Make sure you have the most updated version of P2x  
	* Edit the config file to run the **Ac227TreePlugin**. *See Examples/ADAnalyzer.cfg*   
	* Move to $P2X_ANALYSIS_CODE/ControlScripts   
	* Run `./LaunchBatchScript.py --anacal <Data type> --mode <Config file> --jmax <# of cores>`   
		where <Data type> = WetCommissioning, 180316_Rampdown, or 180316_Background  	
	* The results will be placed in $P2X_ANALYZED   


## Plotting RnPo Events
**(1)** Set all environment variables in your personal .bashrc. *See Examples/ac227examplebashrc.txt*  
	* $P2X_ANALYZED: where the Ac227TreePlugin results live
	* $AD_AC227ANALYSIS_RESULTS: where you want the resulting histograms and TGraphs to live
	* $AD_AC227_PLOTS: where you want the resulting plots to live

**(2)** Create the RNPO class that will contain a TChain of all TAc root trees using the good run list and MakeAcTreeClass.C. 
	* Copy Neutrino18_GoodRuns.txt to the current directory or the most current good run list. For the second case make sure to change the name of the file in the macro.
	* Run the MakeAcTreeClass.C macro. This will read in the list of files from the good run list, find the associated AD1_Extra_Phys.root file in $P2X_ANALYZED, and add the TAc TTree to a TChain.

The next part of the analysis happens in two steps, the calculations and the plotting.
There are four macros that will perform the calculations for the following scenarios: per cell, over time, per column over time, and per row over time.
Each of these macros will create two root files. One that is populated with TGraphs for variables such as rate, Po energy mean, Po energy sigma, etc.
The other is populated with histograms of distributions such as time, energy, psd, etc., for either each cell or each time bin. 

**(3)** Calculate the RnPo rate and other interesting variables, saving information to root trees for future plotting. 
	* Set the desired cuts for PSD, energy [MeV], and position [mm] in RnPoCalculate.sh
	* Set the desired time bin [hrs] in RnPoCalculate.sh
	* Set the desired range start for the dt distribution fit such that x_0 = dtFit x binwidth
	* Run RnPoCalculate.sh



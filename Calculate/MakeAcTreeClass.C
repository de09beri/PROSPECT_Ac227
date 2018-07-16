#include "TChain.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;

int makeAcTreeClass(){
	clock_t tStart = clock();	
	
	ifstream file;
	file.open("/home/dberish/AD_Ac227Analysis/Analysis/Neutrino18_GoodRuns.txt", ifstream::in);
//	file.open("GoodRuns_ReactorOn.txt", ifstream::in);
	
	if(!(file.is_open() && file.good())){
		printf("Good runs file not found. Exiting. \n");
		return -1;
	}


	TChain *chain = new TChain("Ac227TreePlugin/TAc");

	int i=0;	
	while(file.good() & !file.eof()){
		string line;
		getline(file, line);
		TString str = Form("%s/%s/ADAnalyzer.root",gSystem->Getenv("P2X_ANALYZED"),line.data());

		chain->Add(str.Data());
		file.peek();

		if(i%100==0) printf("At file %i: %s \n",i,(const char*)str);
		i++;
	}

	int nadded = chain->GetListOfFiles()->GetEntries();
	chain->Lookup(1);
	int nmissing = nadded - chain->GetListOfFiles()->GetEntries();

	printf("Attempted to add %d files from good run list to Ac227 TChain \n",nadded);
	printf("Unable to add %d files to TChain \n",nmissing);

	printf("Making class \n");
	chain->MakeClass("RNPO");

	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	return 0;
}
